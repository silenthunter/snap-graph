
#include "graph_modularity.h"


void computeModularityValue(graph_t *G, attr_id_t *membership, attr_id_t numCommunities, double *modularity)
{
	attr_id_t i,j;
	attr_id_t comm;
	double mod=0.0;
	double degree_u,degree_v;
	
	for(i=0; i<G->n; i++)
	{
		comm = membership[i];
		degree_u = (double)G->numEdges[i+1] - G->numEdges[i];
		for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
		{
			if(membership[G->endV[j]] == comm)
				mod +=1.0;
		}
		for(j=0; j<G->n; j++)
		{
			degree_v = (double)G->numEdges[j+1] - G->numEdges[j];
			if(comm == membership[j])
				mod -= (double)(degree_u*degree_v)/(double)(2.0*G->m);
		}
	}		
	printf("Modularity is:%g\n",mod/(2*G->m));
}


 void modularity_spectral(graph_t *G, attr_id_t **membership, attr_id_t *numCommunities)
 {
 	attr_id_t *v2C, *degree, *vertex, *v2pos;
	attr_id_t curCommunity=0, newCommunity,toSplit;
	attr_id_t n=G->n,sumV1,sumV2,comm,count1,count2;
	attr_id_t *memblock;
	attr_id_t i,j,communitySize,degreeSum,start;
	list *Q;
	node *first;
	double *eigenVectorOld,*eigenVectorNew;
	

	*numCommunities = 1;

	printf("Allocating memory\n");
	//Allocating memories to the above datastructures
	memblock = (attr_id_t*)malloc(sizeof(attr_id_t)*4*n);
	assert(memblock!=NULL);
	v2C=memblock;			//v2C is a map from vertex to the community it belongs.
	vertex=memblock + 1*n;		//vertex is an array of vertices belonging to a particular community.
	v2pos=memblock + 2*n;		//v2pos is a map from vertex to its respective position in the community.
	degree = memblock + 3*n;	//degree is a map from vertex to its degree in its respective community.
	eigenVectorOld=(double*)malloc(sizeof(double)*n);
	eigenVectorNew=(double*)malloc(sizeof(double)*n);
	assert(eigenVectorOld!=NULL);assert(eigenVectorNew!=NULL);

	printf("Allocating memory done\n");

	for(i=0; i< G->n; i++)
	{
		v2C[i] = 0;
		vertex[i] = i;
		v2pos[i] = i;
		degree[i] = G->numEdges[i+1] - G->numEdges[i];
	}
	//print_attr_id_t_Vector(vertex,0,n);
	*membership=v2C;
	//Making a queue. This queue will store all the communities that are yet to be processed.
	Q=(list*)makeList();
	append(Q, makeNode(curCommunity));
	printList(Q);
	
	//printf("v2C:");print_attr_id_t_Vector(v2C,0,G->n);
	//printf("degree");print_attr_id_t_Vector(degree,0,G->n);
	//printf("vertex:");print_attr_id_t_Vector(vertex,0,G->n);
	//printf("v2pos:");print_attr_id_t_Vector(v2pos,0,G->n);

	while(Q->size > 0)
	//while(Q->size < 4)
	{
		first = (node*) getFirst(Q);
		curCommunity = first->id;
		deleteFirst(Q);

		printf("\n\nEvaluating Community:%d\n",curCommunity);
		communitySize=degreeSum=0;
		//Checking which all vertices belong to this community and updating the vertex Vector accordingly.
		for(i=0; i<G->n; i++)
		{
			if(v2C[i] == curCommunity)
			{
				vertex[communitySize] = i;
				v2pos[i] = communitySize;
				degreeSum += G->numEdges[i+1]-G->numEdges[i];
				communitySize++;
			}
			//eigenVectorNew[i]=pow(-1,i);
		}
		if(communitySize == 1)	continue;

		 computeEigen(G,eigenVectorOld, eigenVectorNew, v2C,v2pos,degree,vertex,&toSplit,curCommunity,communitySize,degreeSum);
		 //printDoubleVector(eigenVectorOld,0,communitySize);

		 //Q->size=-10;

		if(toSplit == 0)
			continue;

		newCommunity = *numCommunities;
		count1=count2=sumV1=sumV2=0;

		for(i=0; i<communitySize; i++)
		{
			if(eigenVectorOld[i] > 0) count1++;
			else count2++;
		}
		if(count1 == 0 || count2 == 0)
			continue;				//All eigen values are of same size and hence no division is required.
		for(i=0; i<communitySize; i++)
		{
			if(eigenVectorOld[i] > 0)
			{
				v2C[vertex[i]] = newCommunity;
				sumV1++;
			}
			else
				sumV2++;
		}
		//printf("sumV1=%d,sumV2=%d,communitySize=%d\n",sumV1,sumV2,communitySize);
		//Now updating the degree Vectors
		for(i=0; i<communitySize; i++)
		{
			comm = v2C[vertex[i]];
			degree[vertex[i]] = 0;
			for(j=G->numEdges[vertex[i]]; j<G->numEdges[vertex[i]+1]; j++)
			{
				if(v2C[G->endV[j]] == comm)
					degree[vertex[i]]++;
			}
		}
		*numCommunities = *numCommunities + 1;
		append(Q,makeNode(curCommunity));
		append(Q,makeNode(newCommunity));

	}
		//printf("v2C:");print_attr_id_t_Vector(v2C,0,G->n);
		//printf("degree:");print_attr_id_t_Vector(degree,0,G->n);
		//printf("vertex:");print_attr_id_t_Vector(vertex,0,G->n);
		//printf("v2pos:");print_attr_id_t_Vector(v2pos,0,G->n);


}



void computeEigen(graph_t *G, double *eigenVectorOld, double *eigenVectorNew, attr_id_t *v2C, attr_id_t *v2pos, int* degree, attr_id_t *vertex, attr_id_t *toSplit, attr_id_t currCommunity, attr_id_t communitySize, attr_id_t degreeSum)
{
	attr_id_t i,j,k;
	attr_id_t iterCount,count=0, niter = communitySize > 1000?communitySize:100;
	niter = 10*log(communitySize);
	double eigenValue,degree_u,degree_v;
	double normalizedSum,ktx,mneg=0.0;

/*	//checking array
	FILE *fp = fopen("arvind.txt", "a+");
	fprintf(fp, "degreeSum=%d\n",degreeSum);
	//printf("v2C:");print_attr_id_t_Vector(v2C,0,G->n);
	//printf("degree:");print_attr_id_t_Vector(degree,0,G->n);
	//printf("vertex:");print_attr_id_t_Vector(vertex,0,communitySize);
	//printf("v2pos:");print_attr_id_t_Vector(v2pos,0,G->n);
	double **arr;
	double value;
	arr = (double**)malloc(sizeof(double*)*communitySize);
	for(i=0;i <communitySize; i++)
	{
		arr[i] = (double*)calloc(communitySize, sizeof(double) );
	}
	
	double  num=0;

	for(i=0; i<communitySize; i++)
	{
		degree_u = G->numEdges[vertex[i]+1] - G->numEdges[vertex[i]];
		num=0;
		value = degree[vertex[i]] - (double)(degree_u*degreeSum)/(double)(2.0*G->m);
		for(j=G->numEdges[vertex[i]]; j<G->numEdges[vertex[i]+1];j++)
		{
			
			if(v2C[G->endV[j]]==currCommunity)
			{
				num++;
				arr[i][v2pos[G->endV[j]]] += (double) (G->dbl_weight_e[j]);
			}

		}
		for(k=0; k< communitySize; k++)
		{
			degree_v = G->numEdges[vertex[k]+1] - G->numEdges[vertex[k]];
			//arr[i][v2pos[vertex[k]]] -=  (double)(degree_u * degree_v)/(double)(2*G->m) + value;
			arr[i][v2pos[vertex[k]]] -=  (double)(degree_u * degree_v)/(double)(2*G->m) ;
		}
		arr[i][i] -=  value;
		
		for(k=0; k<communitySize; k++)
		{
			fprintf(fp, "%g,", arr[i][k]);
		}

		fprintf(fp,"\n");
	}

	fprintf(fp, "\n\n\n\n");
	//checking matrix ends.
*/	while(1)
	{
		iterCount = 0;
		count++;
		normalizedSum=0.0;
		srand(time(NULL));
		for(i=0; i<communitySize; i++)
		{
			eigenVectorOld[i] = 2.0f *((float)rand()/(float)RAND_MAX) -1.0f;
			//eigenVectorOld[i]= 2.0f * (float)drand48() -1.0f;
			normalizedSum += eigenVectorOld[i] * eigenVectorOld[i];
		}
		normalizedSum = sqrt(normalizedSum);
		for(i=0; i<communitySize; i++)
			eigenVectorOld[i] = eigenVectorOld[i]/normalizedSum;

		//printDoubleVector(eigenVectorOld,0,communitySize);

		while(iterCount <niter)
		{
//
			iterCount++;
			ktx=0.0;
			for(i=0; i<communitySize; i++)
			{
				eigenVectorNew[i]=0.0;
				degree_u = G->numEdges[vertex[i]+1] - G->numEdges[vertex[i]];
				for(j=G->numEdges[vertex[i]]; j<G->numEdges[vertex[i]+1]; j++)
				{
					if(v2C[G->endV[j]] == currCommunity)
					{
						eigenVectorNew[i] += G->dbl_weight_e[j] * eigenVectorOld[v2pos[G->endV[j]]];
					}
				}
				eigenVectorNew[i] -= (((double)degree[vertex[i]]) - (double)(degree_u * degreeSum)/(double)(2.0*G->m))* eigenVectorOld[i];
				ktx += (double) degree_u * eigenVectorOld[i];
			}

			ktx /=(double)(2.0 * G->m);
			normalizedSum = 0.0;

			for(i=0; i<communitySize; i++)
			{
				eigenVectorNew[i] -= (double)(G->numEdges[vertex[i]+1]-G->numEdges[vertex[i]]) *ktx;
				eigenVectorNew[i] -= mneg*eigenVectorOld[i];
				normalizedSum += eigenVectorNew[i]*eigenVectorNew[i];
			}
			normalizedSum = sqrt(normalizedSum);

			for(i=0; i<communitySize;i++)
			{
				eigenVectorOld[i] = eigenVectorNew[i]/normalizedSum;
			}
		}
		
		//printf("Old eigenVectoris:");printDoubleVector(eigenVectorOld,0,communitySize);
		eigenValue =0.0;
		ktx=0.0;
		for(i=0; i<communitySize; i++)
		{
			degree_u = G->numEdges[vertex[i]+1] - G->numEdges[vertex[i]];
			for(j=G->numEdges[vertex[i]]; j<G->numEdges[vertex[i]+1]; j++)
			{
				if(v2C[G->endV[j]] == currCommunity)
				{
					eigenValue += G->dbl_weight_e[j] * eigenVectorOld[v2pos[G->endV[j]]] * eigenVectorOld[i];
				}
			}
			eigenValue += -((double)degree[vertex[i]] - (double)(((double) degree_u *degreeSum)/(double)(2.0*G->m)))* eigenVectorOld[i] * eigenVectorOld[i];
			ktx += (double)degree_u * eigenVectorOld[i];
		}

		ktx /=(double)(2.0 * G->m);
		for(i=0; i<communitySize; i++)
			eigenValue -= ((double)(G->numEdges[vertex[i]+1]- G->numEdges[vertex[i]])) *ktx *eigenVectorOld[i];
	
		printf("The eigenValue is %f\n",eigenValue);
		//fprintf(fp, "The eigenValue is %f\n",eigenValue);

		if(eigenValue <0)
		{
			if(count==2)
			{
				*toSplit = 0;
				break;
			}
			mneg=eigenValue;
		}
		else
		{	
			*toSplit = 1;
			break;
		}
	}
	//fclose(fp);
//

/*
		//temporary code::::::::::::::::;
			iterCount++;
			normalizedSum = 0.0;
			for (i=0; i<communitySize; i++)
			{
				eigenVectorNew[i]=0.0;
				for(j=0; j<communitySize; j++)
				{
					eigenVectorNew[i] += eigenVectorOld[j] * arr[i][j];
				}
				eigenVectorNew[i] -= mneg * eigenVectorOld[i];
				normalizedSum += eigenVectorNew[i] * eigenVectorNew[i];
			}
			normalizedSum = sqrt(normalizedSum);
			for(i=0; i<communitySize; i++)
			{
				eigenVectorOld[i] = eigenVectorNew[i] / normalizedSum;
			}
		}	
		printf("Old eigenVectoris:");printDoubleVector(eigenVectorOld,0,communitySize);
		
		eigenValue = 0.0;
		for(i=0; i<communitySize; i++)
		{
			for(j=0; j<communitySize; j++)
			{
				eigenValue += arr[i][j] * eigenVectorOld[j] * eigenVectorOld[i] ;
			}

		}
		printf("The eigenValue is %f\n",eigenValue);
		fprintf(fp, "The eigenValue is %f\n",eigenValue);

		if(eigenValue <0)
		{
			if(count==2)
			{
				*toSplit = 0;
				break;
			}
			mneg=eigenValue;
		}
		else
		{	
			*toSplit = 1;
			break;
		}
	}
	fclose(fp);
*/

}


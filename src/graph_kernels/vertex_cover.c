

#include "graph_vertex_cover.h"




void calculateVertexCover(graph_t *G)
{	

	printf("Hello World123\n");
	
	double *wp_v;					//weight associated with each vertex.
	double *delta_e;				//Delta of each edge as mentioned in the algorithm
	attr_id_t *degree_v;				//Degree of each vertex
	attr_id_t *visited_e, *visited_v;		//Whether that edge has been visited or not. Ditto for vertex. Visited_v is the final cover.
	attr_id_t *position_e;				//Position stores the corresponding position of the undirected edge for each edge.
	attr_id_t i,j,u,v,k, edge_counter;

	double *memblock;
	attr_id_t *memblock1;
	memblock = (double*) malloc(sizeof(double)*(G->n+2*G->m));
	wp_v = memblock;
	delta_e = memblock + G->n;
	memblock1 = (attr_id_t*) malloc(sizeof(attr_id_t)*(2*G->n + 4*G->m));
	degree_v = memblock1;
	visited_v = memblock1 + G->n;
	visited_e = memblock1 + 2*G->n;
	position_e = memblock1 + 2*(G->n + G->m);

	#pragma omp parallel for private(u,v,j,k)
	for(i=0; i<G->n; i++)
	{
		wp_v[i] = G->dbl_weight_v[i];
		degree_v[i] = G->numEdges[i+1] - G->numEdges[i];
		//printf("degree_v= %d\n",degree_v[i]);
		visited_v[i] = 0;
		if(degree_v[i] == 0)
			visited_v[i]=1;
		for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
		{
			u = i;
			v = G->endV[j];
			delta_e[j] = 0;
			visited_e[j] = 0;
			//printf("u=%d, v=%d\n",u,v);
			if(v < u )
				continue;		// we have already covered this case when we visited v.
			for (k=G->numEdges[v]; k<G->numEdges[v+1]; k++)
			{
				if(G->endV[k] == u)
					break;
			}
			position_e[j] = k;
			position_e[k] = j;
		}
	}
	printf("degree:");print_attr_id_t_Vector(degree_v,0,G->n);
	//printf("position:\n");	print_attr_id_t_Vector(position_e,0,2*G->m);

	printf("Starting real calculation\n");	
	edge_counter = 2*G->m;
	double val1,val2;
	int count =0;
	while(edge_counter > 0)
	//while(count <3)
	{
		count ++;
		//printf("%d\t",edge_counter);
		#pragma omp parallel
		{
		#pragma omp for private(j,u,v,val1,val2)
		for (i=0; i< G->n; i++)
		{
		//	printf("count=%d\n",count);
			if (visited_v[i] == 1)
				continue;
			for(j=G->numEdges[i]; j< G->numEdges[i+1]; j++)
			{
				if(visited_e[j] == 1 )
					continue;
				u = i;
				v = G->endV[j];
				//printf("aa%d\t",degree_v[v]);
				val1 = wp_v[u]/degree_v[u];
				val2 = wp_v[v]/degree_v[v];
				delta_e[j] = val1 < val2 ? val1 : val2;
			}
		}
		#pragma omp for private(j) reduction(-:edge_counter)
		for(i=0; i<G->n; i++)
		{
			//printf("This vertex %d is executed by thread %d\n",i,omp_get_thread_num());
			if (visited_v[i] == 1)
				continue;
			double sum = 0.0;
			for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
			{
				if(visited_e[j] == 1)
					continue;
				sum += delta_e[j];
			}
			wp_v[i] -= sum;
			if(wp_v[i] <= 0.00001)	//aka this vertex is in VC.
			{
				visited_v[i] = 1;
				edge_counter -= degree_v[i]*2;	//It is multiplied by because it is an undirected graph.
				for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
				{
					if(visited_e[j] == 1)
						continue;
					visited_e[j] = 1;
					degree_v[G->endV[j]] -= 1;
					visited_e[position_e[j]] = 1;
				}
				
			}
		}
		}
		
	}
	print_attr_id_t_Vector(visited_v,0,G->n);
	count=0;
	for(i=0; i<G->n; i++)
	{
		if(visited_v[i] == 1 && ((G->numEdges[i+1] - G->numEdges[i])>0))
			count ++;
	}
	printf("The size of VC = %d\n", count);
}




void calculateUnweightedVertexCover(graph_t *G)
{
	attr_id_t i,j,u,v;
	attr_id_t max, max_e, max_u, max_v,edge_counter;
	
	attr_id_t *visited_v, *visited_e;
	attr_id_t * degree_v;
	attr_id_t *memblock;
	memblock = (attr_id_t*) malloc(sizeof(attr_id_t)*(2*G->m + 2*G->n));
	//memblock = (attr_id_t*) malloc(sizeof(attr_id_t)*(G->m + 2*G->n));
	visited_v = memblock;
	degree_v = memblock + G->n;
	visited_e = memblock + 2*G->n;
	
	#pragma omp parallel 
	{
		#pragma omp for
		for(i=0; i<G->n; i++)
		{
			visited_v[i] = 0;
			degree_v[i] = G->numEdges[i+1] - G->numEdges[i];
		}
		#pragma omp for
		//for(i=0; i<G->m; i++)
		for(i=0; i<2*G->m; i++)
		{
			visited_e[i] = 0;
		}
	}
	//print_attr_id_t_Vector(degree_v,0,G->n);
	edge_counter = 2*G->m;
	//edge_counter = G->m;
	while(edge_counter > 0)
	{
		max = 0;
		#pragma omp parallel for shared(max,max_e,max_u,max_v) private(j,u,v)
		for(i=0; i<G->n; i++)
		{
			if(degree_v[i] == 0)
				continue;
			for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
			{
				u = i;
				v = G->endV[j];
				if(degree_v[u] + degree_v[v] > max)
				{
					max = degree_v[u]+ degree_v[v];
					max_e = j;
					max_u = u;
					max_v = v;
				}
			}
		}
		edge_counter -= max;
		//printf("edgecounter=%d\n",edge_counter);
		visited_e[max_e] =1;
		degree_v[max_u] = 0;
		visited_v[max_u] = 1;
		degree_v[max_v] = 0;
		visited_v[max_v] = 1;
		//print_attr_id_t_Vector(visited_e,0,2*G->m);
		//print_attr_id_t_Vector(visited_v,0,G->n);

	}
	attr_id_t count = 0;
	for(i=0; i<G->n; i++)
	{
		if(visited_v[i] == 1)
			count ++;
	}
	printf("The size of VC = %d\n", count);
}

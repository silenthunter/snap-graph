#include "graph_defs.h"
#include "graph_gen.h"
#include "network.h"


void network_to_graph(graph_t *G, NETWORK *N);

void read_GML_graph(graph_t* G, char* filename) 
{

	NETWORK *N = (NETWORK*)malloc(sizeof(NETWORK));;
	FILE *ifile = fopen(filename, "r");
	printf("Reading from file...\n");
	read_network(N,ifile);	//Note that readgml is hardcoded to read undirected graphs,even if the input graph is undirected.
	printf("Reading from file done\n");
	network_to_graph(G,N);
}	

void network_to_graph(graph_t *G, NETWORK *N)
{
	G->n = N->nvertices;
	int i,numEdges=0;
	VERTEX vertex;
	int j,degree,target,sumDegree=0;
	double weight;
	attr_id_t start;
	for(i=0; i<N->nvertices; i++)
	{
		vertex = N->vertex[i];
		numEdges+=vertex.degree;
	}
	assert(numEdges%2==0);
	G->m = numEdges/2;
	assert(G->n > 0);
	assert(G->m > 0);
	G->undirected=1;	//Hard-coded undiected
	G->zero_indexed=0;	//Hard-coded false
	G->weight_type=4;	//Hard-coded weight_type is double
	//Allocating memory 
	G->numEdges = (attr_id_t*) calloc(G->n+1, sizeof(attr_id_t) );
        G->endV = (attr_id_t*) calloc(2*G->m, sizeof(attr_id_t)  );
        G->dbl_weight_e = (double*) malloc(sizeof(double)* 2*G->m  );
        if(G->numEdges==NULL || G->endV==NULL || G->dbl_weight_e==NULL)
        {
                printf("Error in memory alloc: Requested memory size = %d\n",sizeof(attr_id_t)*(2*G->m+ G->n) + sizeof(double)*G->m*2 );
                exit;
        }
	G->numEdges[0]=0;
	int count=0;
	for(i=0; i<N->nvertices; i++)
	{
		vertex = N->vertex[i];
		degree = vertex.degree;
		G->numEdges[i+1] = G->numEdges[i] + degree;
		start = G->numEdges[i];

		for(j=0; j<degree;j++)
		{
			target = vertex.edge[j].target;
			weight = vertex.edge[j].weight;
			G->endV[start+j] = target;
			G->dbl_weight_e[start+j] = weight;
			count++;
		}
	}
	free_network(N);
	free(N);

} 

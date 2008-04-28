

#include "graph_vertex_cover.h"

void readSample(graph_t*);

int main(int argc, char** argv)
{
	printf("Starting vertex cover\n");
	graph_t G;
	attr_id_t i;
	//readSample(&G);

	//read_DIMACS_graph(&G,"vcdata/frb59-26-mis/frb59-26-2.mis");
	read_DIMACS_graph(&G,"vcdata/frb56-25-mis/frb56-25-1.mis");
	//read_DIMACS_graph(&G,"vcdata/frb100-40.mis");
	// gen_RMAT_graph(&G, "./test.rmat");
	G.dbl_weight_v = (double*)malloc(sizeof(double)*G.n);
	
	G.m /=2;

	for(i=0; i<G.n; i++)
	{
		
		G.dbl_weight_v[i] = 1;
	}
	//print_graph(&G);
	printf("reading graph done\n");
	double startTime, endTime;
	startTime = omp_get_wtime();
	calculateVertexCover(&G);
	endTime = omp_get_wtime();
	printf("Time taken %g\n",endTime - startTime);
	
	startTime = omp_get_wtime();
	//calculateUnweightedVertexCover(&G);
	endTime = omp_get_wtime();
	printf("Time taken %g\n",endTime - startTime);
}






attr_id_t endV[] = {1,2,0,2,0,1,3,2,4,5,3,6,3,6,4,5};
attr_id_t numEdges[]= {0,2,4,7,10,12,14,16};
double dbl_weight_e[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
double dbl_weight_v[] = {1,2,1,2,1,2.0,1.0};

void readSample(graph_t *G) 
{
        G->n = 7;
        G->m = 8;
	G->endV = endV;
	G->numEdges = numEdges;
	G->weight_type=4;
	G->zero_indexed=1;
	G->dbl_weight_e = dbl_weight_e;
	G->dbl_weight_v = dbl_weight_v;
	G->undirected=1;
}



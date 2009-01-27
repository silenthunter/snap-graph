#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "graph_partitioning.h"
#include "utils.h"

int main()
{
	printf("\n Starting driver of modularity\n");
	graph_t G;
	int *membership, numCommunities,i;
	double modularity, startTime, endTime;
	
	read_GML_graph(&G,"./gml/karate.gml");
	//gen_RMAT_graph(&G, "./test.rmat"); G.m/=2;	
	//readSample(&G);
    //	read_DIMACS_graph(&G, "./gml/celegans_metabolic_modified.net");G.m/=2;
    //	read_DIMACS_graph(&G, "./gml/email_modified.txt");G.m/=2;
    //	read_DIMACS_graph(&G, "./gml/jazz_modified.net");G.m/=2;
    //	read_DIMACS_graph(&G, "./gml/www.dat");G.m/=2;
	//G.dbl_weight_e = (double*)malloc(sizeof(double)*2*G.m);
	//for(i=0; i<2*G.m;i++)
	//	G.dbl_weight_e[i]=1.0;
	printf("reading graph done%d\n",G.n);
    //	print_graph(&G);

	startTime = get_seconds();
	//modularity_spectral(&G,&membership,&numCommunities,1);
	modularity_spectral_wo_klin(&G,&membership,&numCommunities);
	//modularity_spectral(&G,&membership,&numCommunities);
	endTime = get_seconds();
	//int comms[] = {1,1,1,1,3,3,3,1,0,0,3,3,1,1,0,0,3,1,0,1,0,1,0,2,2,2,0,2,2,0,0,2,0,0};
	int comms[]={4,5,5,5,2,2,2,5,1,5,2,2,4,5,1,1,2,4,1,4,1,4,1,3,3,3,1,3,3,1,1,3,1,1};
	printf("Time taken = %g\n",endTime - startTime);
	
	computeModularityValue(&G,membership,numCommunities,&modularity);
	//computeModularityValue(&G,comms3,4,&modularity);
	printf("Num communities=%d\n",numCommunities);

}

attr_id_t endV[] = {1,2,0,2,0,1,3,2,4,5,3,6,3,6,4,5};
attr_id_t numEdges[]= {0,2,4,7,10,12,14,16};
double dbl_weight_e[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

void readSample(graph_t *G)
{
	G->n = 7;
	G->m = 8;
	G->endV = endV;
	G->numEdges = numEdges;
	G->weight_type=4;
	G->zero_indexed=1;
	G->dbl_weight_e = dbl_weight_e;
	G->undirected=1;
}

	





#include "graph_defs.h"
#include "graph_metrics.h"
#include "utils.h"
#include "sprng.h"
#include "bc_cuda.h"

void cuda_betweenness_centrality(graph_t* G, double* BC2, long numSrcs){
	int numVert;
	int numEdge;
	int* edges;
	float* BC;
	int* glob;
	float* dep;
	cuda_bc(numVert, numEdge, edges, BC, glob, dep);
}

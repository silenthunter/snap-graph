#include "utils.h"

double get_seconds() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return (double) (tp.tv_sec + ((1e-6)*tp.tv_usec));
}

void usage_graph_options() {
    fprintf(stdout, "The input file name is a required argument. We try to determine the graph");
    fprintf(stdout, " format from the file extension; optionally, you can specify the graph type with");
    fprintf(stdout, " the -graph option.\n\n");
    fprintf(stdout, "Supported graph file formats (for use with the -graph option):\n");
    fprintf(stdout, "snap    (.gr)       SNAP file format (default file format supported by this package).\n");
    fprintf(stdout, "dimacs  (.dim)      graph format used in the 9th DIMACS Shortest Paths Challenge.\n");
    fprintf(stdout, "metis   (.met)      graph representation used by the Metis partitioning package.\n");
    fprintf(stdout, "graphml (.graphml)  GraphML representation.\n");
    fprintf(stdout, "gml     (.gml)      GML format.\n");
    fprintf(stdout, "rand    (.rnd)      a synthetic random graph generator is invoked."
            " The input file specifies the generator parameters.\n");
    fprintf(stdout, "rmat    (.rmat)       a synthetic scale-free graph generator is invoked."
            " The input file specifies the generator parameters.\n");
    fprintf(stdout, "sqmesh  (.sqm)      a synthetic 2D square mesh generator is invoked."
            " The input file specifies the generator parameters.\n");
    fprintf(stdout, "lmesh   (.lm)       a synthetic 2D long mesh generator is invoked."
            " The input file specifies the generator parameters.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "The default output filename is \"results.out\"."
            " Specify an alternate file name with the -outfile option.\n");
}

void graph_ext_check(char* filename, char* graph_type) {

    char *ptr, *tmp; 

    ptr = filename;
    
    /* Find the last "." in the file name */
    do {
        tmp = strstr(ptr+1, ".");
        if (tmp != NULL) 
            ptr = tmp;
    }  while (tmp != NULL);
    ptr++;   

    if (graph_type[0] == 0) {
        if (strcmp(ptr, "gr") == 0) {
            strcpy(graph_type, "snap");
        }

        if (strcmp(ptr, "dim") == 0) {
            strcpy(graph_type, "dimacs");    
        }

        if (strcmp(ptr, "met") == 0) {
            strcpy(graph_type, "metis");    
        }

        if (strcmp(ptr, "graphml") == 0) {
            strcpy(graph_type, "graphml");    
        }

        if (strcmp(ptr, "gml") == 0) {
            strcpy(graph_type, "gml");    
        }

        if (strcmp(ptr, "rand") == 0) {
            strcpy(graph_type, "rand");    
        }

        if (strcmp(ptr, "rmat") == 0) {
            strcpy(graph_type, "rmat");    
        }

        if (strcmp(ptr, "sqmesh") == 0) {
            strcpy(graph_type, "sqmesh");    
        }

        if (strcmp(ptr, "lmesh") == 0) {
            strcpy(graph_type, "lmesh");    
        }
    }

}

void prefix_sums(attr_id_t *input, attr_id_t* result, attr_id_t* p, long n) {

#ifdef _OPENMP
    attr_id_t i, j, r, start, end, add_value;
    int tid, nthreads;

    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    r = n/nthreads;

    result[0] = 0;

    #pragma omp for
    for (i=1; i<n+1; i++)
        result[i] = input[i-1];

    start =  tid*r + 1;
    end   = (tid+1)*r;

    if (tid == nthreads-1)
        end = n+1;

    for (j=start; j<end; j++)
        result[j] = input[j-1] + result[j-1];

    p[tid] = result[end-1];

    #pragma omp barrier

    if (tid == 0) {
        for (j=1; j<nthreads; j++)
            p[j] += p[j-1];
    }

    #pragma omp barrier

    if (tid>0) {
        add_value=p[tid-1];
        for (j=start-1; j<end; j++)
            result[j] += add_value;
    }

    #pragma omp barrier
#else

    result[0] = 0;
    for (int i=1; i<n+1; i++) {
        result[i] = result[i-1] + input[i-1];
    }
#endif
    
}



void print_graph(graph_t* G)
{
	int i,j;
	int degree,start;
	
	printf("\n Printing graph\n");
	printf("Number of vertices =%d, edges=%d\n",G->n,G->m);
	if(G->undirected) 
		printf("Graph is undirected\n");
	else
		printf("Graph is directed\n");
	printf("Graph weighttype is %d\n",G->weight_type);


	for(i=0;i<G->n;i++)
	{
		degree=G->numEdges[i+1]-G->numEdges[i];
		start=G->numEdges[i];
		printf("Degree of vertex i=%d is %d. Its neighbours are:",i,degree);
		
		for(j=0;j<degree;j++)
		{
			if(G->weight_type == 1) 
				printf("%d[%d], ",G->endV[start+j], G->int_weight_e[start+j]);
			else if(G->weight_type == 2) 
				printf("%d[%d], ",G->endV[start+j], G->l_weight_e[start+j]);
			else if(G->weight_type == 3) 
				printf("%d[%f], ",G->endV[start+j], G->fl_weight_e[start+j]);
			else if(G->weight_type == 4) 
				printf("%d[%f], ",G->endV[start+j], G->dbl_weight_e[start+j]);
		}
		printf("\n");
	}
}


	





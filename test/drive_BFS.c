#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "utils.h"

int main(int argc, char** argv) {

    char *infilename, *outfilename, *graph_type;
    FILE* fp;    
    graph_t* g;

    long src;
    int curArgIndex;
    long numSrcs;
    int est_diameter;
    long i, j;

    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s (-src <vertex ID (0 to n-1)>) -infile <graph filename>"
                " (-graph <graph type> -outfile <output filename>)\n\n", "eval_BFS");
        
        usage_graph_options();
        
        fprintf(stdout, "Using the -src option, specify the vertex ID to run BFS from. A random vertex is selected if the src is not specified.\n\n");
        exit(-1);
    }

    curArgIndex = 0;
    infilename = (char *) calloc(500, sizeof(char));
    outfilename = (char *) calloc(500, sizeof(char));
    graph_type = (char *) calloc(500, sizeof(char));

    strcpy(outfilename, "/tmp/results.out");

    src = -1;

    while (curArgIndex < argc) {
        
        if (strcmp(argv[curArgIndex],"-src")==0) {
            src = atol(argv[++curArgIndex]);
        }

        if (strcmp(argv[curArgIndex],"-infile")==0) {
            strcpy(infilename, argv[++curArgIndex]);
        }

        if (strcmp(argv[curArgIndex], "-outfile")==0) {
            strcpy(outfilename, argv[++curArgIndex]);
        } 
 
        if (strcmp(argv[curArgIndex], "-graph")==0) {
            strcpy(graph_type, argv[++curArgIndex]);
        } 
        curArgIndex++; 
    }

    fp = fopen(infilename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error! Could not open input file. Exiting ...\n");
        exit(-1);
    }
    fclose(fp);

    fp = fopen(outfilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error! Could not write to output file. Exiting ...\n");   
        exit(-1);
    }
    fclose(fp);

    graph_ext_check(infilename, graph_type);

    fprintf(stdout, "\n");
    fprintf(stdout, "Input Graph File    : %s\n", infilename);
    fprintf(stdout, "Output Graph File   : %s\n\n", outfilename);



    /* Step 2: Generate graph */
    g = (graph_t *) malloc(sizeof(graph_t));
    graph_gen(g, infilename, graph_type);
   
    fprintf(stdout, "No. of vertices     : %ld\n", g->n);
    fprintf(stdout, "No. of edges        : %ld\n\n", g->m);
   
    if (src == -1)
       src = lrand48() % g->n;

    assert((src >= 0) && (src < g->n));
    fprintf(stdout, "Source vertex       : %ld\n\n", src);

    /* Step 3: Run algorithm */

    /* Assuming a low diameter graph */
    est_diameter = 100;
    BFS_parallel_frontier_expansion(g, src, est_diameter);
 
    /* Step 4: Clean up */
    free(infilename);
    free(outfilename);
    free(graph_type);

    free_graph(g);
    free(g);

    return 0;
}

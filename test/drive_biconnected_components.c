#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "utils.h"

int main(int argc, char** argv) {

    char *infilename, *outfilename, *graph_type;
    FILE* fp;    
    graph_t* g;

    int curArgIndex;
    long i, j;
    
    attr_id_t* component_num;

    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s -infile <graph filename>"
                " (-graph <graph type> -outfile <output filename>)\n\n", "eval_biconnected_components");
        
        usage_graph_options();
        exit(-1);
    }

    curArgIndex = 0;
    infilename = (char *) calloc(500, sizeof(char));
    outfilename = (char *) calloc(500, sizeof(char));
    graph_type = (char *) calloc(500, sizeof(char));

    strcpy(outfilename, "/tmp/results.out");

    while (curArgIndex < argc) {

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

    /* Step 3: Run algorithm */ 

    component_num = (attr_id_t *) malloc(g->m*sizeof(attr_id_t));
    biconnected_components(g, component_num);
 
    /* Step 4: Clean up */
    free(component_num);

    free(infilename);
    free(outfilename);
    free(graph_type);

    free_graph(g);
    free(g);

    return 0;
}

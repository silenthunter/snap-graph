#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "graph_partitioning.h"
#include "utils.h"

attr_id_t* getSeeds(attr_id_t *n_seeds, char* seedsfile) {
	attr_id_t* seeds;
	FILE* fp;
	attr_id_t seed, i;

	fp = fopen(seedsfile, "r");
	fscanf(fp, "%ld", n_seeds);
	if(n_seeds <=0)
		return;

	seeds = (attr_id_t*) calloc (*n_seeds, sizeof(attr_id_t));
	i = 0;
	while(fscanf(fp, "%ld", &seed) > 0) {
		if(i >= *n_seeds) {
			fprintf (stderr, "Error! No. of seeds in the file is"
					 " greater than the specified number %d. "
					 " Exiting ...\n", *n_seeds);
			exit(-1);
		}
		seeds[i++] = seed;
	}
	if (i < *n_seeds) {
		fprintf (stderr, "Error! No. of seeds in the file is"
				 " less than the specified number %ld. Exiting ...\n", *n_seeds);
		exit(-1);
	}
	
	fclose(fp);
	return seeds;
}

int main(int argc, char** argv) {
    char *infilename, *outfilename, *seedsfile, *graph_type, *alg_type;
    FILE* fp;
    graph_t* g;
    int curArgIndex;
	int opt_seed, opt_seeds;
    attr_id_t *membership, *seeds;
	attr_id_t seed;
	attr_id_t n_seeds;
	attr_id_t commId;
    attr_id_t num_visited;
    attr_id_t ext_edges, int_edges;
	int uniform;
    int num_communities;
    double modularity, mod_val, time0;
    long i, size;
    int run_approxBC;
    int digits;

    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s -infile <graph filename>"
                " -graph <graph type> -outfile <output filename>)"
                " -alg <algorithm type>"
				" (-seed <seed vertex> | -seeds <seeds filename>"
				" -uniform <0|1>)"
				"\n\n",
                "eval_seed_community_detection");
        fprintf(stdout, "Algorithm type can be one of the following:\n"
                "CNM         -- greedy agglomerative strategy of "
                "Clauset, Newman and Moore.\n"
                "MB          -- McCloskey and Bader: "
                "normalize CNM Qij's by stddev(Qij).\n"
                "RAT         -- normalize CNM by size ratios.\n"
                "MBRAT       -- normalize McCloskey-Bader by size ratios.\n"
                "LIN         -- McCloskey and Bader: weighted graphs "
                "with linear model.\n"
#if 0
                "WT1/WT2/WT3 -- Wakita and Tsurami's"
                "consolidation ratio heuristics.\n"
                "DDA         -- Normalized modularity heuristics "
                "by Danon, Diaz-Guilera, and Arenas.\n"
                "CCA         -- An agglomerative clustering heuristic "
                "based on local clustering coefficient.\n"
#endif
                "\n");

        usage_graph_options();
        exit(-1);
    }

	opt_seed = 0;
	opt_seeds = 0;
	seed = -1;
	n_seeds = 0;
	uniform = 0;
	commId = -1;
    curArgIndex = 0;
    run_approxBC = 0;
    infilename = (char *) calloc(500, sizeof(char));
    outfilename = (char *) calloc(500, sizeof(char));
	seedsfile = (char *) calloc(500, sizeof(char));
    graph_type = (char *) calloc(500, sizeof(char));
    alg_type = (char *) calloc(500, sizeof(char));
	seeds = NULL;
    size  = 0;

    strcpy(alg_type,"CNM");
    strcpy(outfilename, "output.txt");

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

        if (strcmp(argv[curArgIndex],"-alg")==0) {
#ifdef DEBUGXX
            printf("argv: %s (len: %d)\n",argv[curArgIndex+1], strlen(argv[curArgIndex+1]));
#endif
            strcpy(alg_type, argv[++curArgIndex]);
            printf("alg_type: %s\n",alg_type);
        }
		
		if (strcmp(argv[curArgIndex], "-seed")==0) {
			seed = atoi(argv[++curArgIndex]);
			opt_seed = 1;
		}
		
		if (strcmp(argv[curArgIndex], "-seeds")==0) {
			strcpy(seedsfile, argv[++curArgIndex]);
			opt_seeds = 1;
		}
		
		if (strcmp(argv[curArgIndex], "-uniform")==0) {
			uniform = atoi (argv[++curArgIndex]);
			uniform = (uniform > 0)? 1 : 0;
		}

        curArgIndex++; 
    }

    fp = fopen(infilename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error! Could not open input file. Exiting ...\n");
        exit(-1);
    }
    fclose(fp);
	
	if (opt_seeds == 1) {
		if (opt_seed == 1) {
			fprintf(stderr, "Warning! Both -seed and -seeds specified.\n");
			fprintf(stderr, "Warning! Ignoring -seed option.\n");
			opt_seed = 0;
		}
		fp = fopen(seedsfile, "r");
		if (fp == NULL) {
			fprintf(stderr, "Error! Could not open seed input file. Exiting ...\n");
			exit(-1);
		}
		fclose(fp);
	}
	else if(opt_seed == 0) {
		fprintf (stderr, "Error! -seed or -seeds option is required.\n");
		exit(-1);
	}

    fp = fopen(outfilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error! Could not write to output file. Exiting ...\n");   
        exit(-1);
    }
    fclose(fp);

    graph_ext_check(infilename, graph_type);

    fprintf(stdout, "\n");
    fprintf(stdout, "Input Graph File    : %s\n", infilename);
	if(opt_seeds == 1) {
		fprintf(stdout, "Seeds Input File    : %s\n", seedsfile);
	}
    fprintf(stdout, "Output Graph File   : %s\n\n", outfilename);
	fprintf(stdout, "Algorithm type         : %s\n", alg_type);

	/* Step 2: Read seed input */
	if(opt_seeds == 1) {
		seeds = getSeeds(&n_seeds, seedsfile);
		fprintf (stdout, "Number of seeds        : %ld\n", n_seeds);
	}

    /* Step 3: Generate graph */
    g = (graph_t *) malloc(sizeof(graph_t));
    graph_gen(g, infilename, graph_type);

    fprintf(stdout, "Number of vertices     : %ld\n", g->n);
    if (g->undirected)
        fprintf(stdout, "Number of edges        : %ld\n\n", g->m/2);
    else 
        fprintf(stdout, "Number of edges        : %ld\n\n", g->m);

    if (g->undirected == 0) {
        fprintf(stderr, "Error: the graph has to be undirected.\n");
        fprintf(stderr, "Please check input file.\n");
        exit(-1);
    }

    /* Step 4: Run algorithm */

    membership = (attr_id_t *) malloc(g->n * sizeof(attr_id_t));
    assert(membership != NULL);
    modularity = 0.0;
    num_communities = 0;
    time0 = get_seconds();
	if (opt_seeds == 1) {
		if (uniform)
			fprintf (stderr, "Adopting uniform strategy (arg 'uniform = 1')\n");
		else
			fprintf (stderr, "Adopting greedy strategy (arg 'uniform = 0')\n");

	    seed_set_community_detection(g, alg_type, seeds, n_seeds, uniform,
				membership, &num_communities, &modularity);
	}
	else {
		single_seed_community_detection (g, alg_type, seed, membership,
				&num_communities, &modularity);
	}
    time0 = get_seconds() - time0;

    mod_val = get_community_modularity(g, membership, num_communities);
	
    if (opt_seeds == 1) {
		commId = membership[seeds[0]];
	}
	else {
		commId = membership[seed];
	}
    int_edges = community_edges_internal (g, membership, commId);
    ext_edges = community_edges_external (g, membership, commId);

    /* Step 4: Write output to file */
    fp = fopen(outfilename, "w");
    fprintf(fp, "Number of communities: %d\n", num_communities);
    fprintf(fp, "Modularity score: %f (full), %f (w/o dup)\n", mod_val,
            modularity);
    fprintf(fp, "\n<Vertex ID> <Community ID>\n\n");

	for (i=0; i<g->n; i++) {
		if (commId == membership[i]) {
            size++;
			if (g->zero_indexed)
				fprintf (fp, "%ld\t%d\n", i, membership[i]);
			else
				fprintf (fp, "%ld\t%d\n", i, membership[i]);
		}
	}

    fprintf(stderr, "Number of communities: %d\n", num_communities);
    fprintf(stderr, "Modularity: %f (full), %f (w/o dup)\n", mod_val,
            modularity);
    fprintf(stderr, "Number of internal edges: %ld\n", int_edges);
    fprintf(stderr, "Number of external edges: %ld\n", ext_edges);
    fprintf(stderr, "Execution Time: %f seconds\n", time0);

    /* BFS community */
    if (opt_seeds != 1) {
        fprintf (stderr, "\n*** BFS community detection. ***\n");

        time0 = get_seconds() - time0;
        num_visited = BFS_Community(g, seed, size, membership, &num_communities);
        time0 = get_seconds() - time0;

        mod_val = get_community_modularity(g, membership, num_communities);
        int_edges = community_edges_internal (g, membership, commId);
        ext_edges = community_edges_external (g, membership, commId);

        fprintf(fp, "\n\nBreadth first community :-");
        fprintf(fp, "\nNumber of communities: %d", num_communities);
        fprintf(fp, "\nModularity score: %f (w/o dup)\n", mod_val);
        fprintf(fp, "\n<Vertex ID> <Community ID>\n\n");
		commId = membership[seed];
        for (i = 0; i < g->n; i++) {
		    if (commId == membership[i]) {
    			if (g->zero_indexed)
				    fprintf (fp, "%ld\t%d\n", i, membership[i]);
			    else
			    	fprintf (fp, "%ld\t%d\n", i, membership[i]);
		    }
        }
        
        fprintf(stderr, "Number of communities: %d\n", num_communities);
        fprintf(stderr, "Modularity score: %f (w/o dup)\n", mod_val);
        fprintf(stderr, "Number of internal edges: %ld\n", int_edges);
        fprintf(stderr, "Number of external edges: %ld\n", ext_edges);
        fprintf(stderr, "Execution Time: %f seconds\n", time0);
    }

#if 0
    {
        int *marked;
        long i, j, v;
        long edgecount, eIdx, idx, currLabelCount, currLabelSources;

        typedef struct {
            int from;
            int from_label;
            int to;
            int to_label;
        } cross_edge_t;

        cross_edge_t *edges;

        int cross_edge_compar(const void *v1, const void *v2) {
            const cross_edge_t *p1, *p2;
            p1 = (cross_edge_t *)v1;
            p2 = (cross_edge_t *)v2;

            if ((p1->from == p2->from) && (p1->from_label == p2->from_label) && 
                    (p1->to   == p2->to)   && (p1->to_label   == p2->to_label))
                return 0;

            if (p1->from_label < p2->from_label)
                return -1;

            if (p1->from_label > p2->from_label)
                return +1;

            return (p1->from < p2->from ? -1 : +1);
        }


        marked = (int *)malloc(g->n * sizeof(int));
        assert(marked != NULL);

        for (i=0; i<g->n; i++)
            marked[i] = (g->numEdges[i+1] > g->numEdges[i] ? 1:  0); /* Degree > 0 */

        edgecount = 0;
        for (i=0; i<g->n; i++) {
            if (marked[i]) {
                for (j=g->numEdges[i] ; j<g->numEdges[i+1] ; j++) {
                    v = g->endV[j];
                    if (!marked[v]) fprintf(stderr,"ERROR: v %d not marked\n",v);
                    if (membership[i] != membership[v]) {
#if 0
                        fprintf(stderr,"TEST %4d (%4d) != %4d (%4d)\n",
                                i, membership[i], v, membership[v]);
#endif
                        edgecount++;
                    }
                }
            }
        }

        if (edgecount > 0) {

            edges = (cross_edge_t *)malloc(edgecount * sizeof(cross_edge_t));
            assert(edges != NULL);


            eIdx = 0;
            for (i=0; i<g->n; i++) {
                if (marked[i]) {
                    for (j=g->numEdges[i] ; j<g->numEdges[i+1] ; j++) {
                        v = g->endV[j];
                        if (membership[i] != membership[v]) {
                            edges[eIdx].from = i;
                            edges[eIdx].from_label = membership[i];
                            edges[eIdx].to   = v;
                            edges[eIdx].to_label   = membership[v];
                            eIdx++;
                        }
                    }
                }
            }

            qsort(edges, edgecount, sizeof(cross_edge_t), cross_edge_compar);

#if 1
            for (i=0 ; i<edgecount ; i++)
                fprintf(stderr,"crossedge[%4d]: %4d %4d %4d %4d\n",
                        i, edges[i].from, edges[i].from_label, edges[i].to, edges[i].to_label);
#endif

            fprintf(stdout,"\n");

            idx = 0;
            currLabelCount   = 1;
            currLabelSources = 1;
            for (i=1 ; i<edgecount ; i++) {
                if (edges[i].from_label == edges[idx].from_label) {
                    currLabelCount++;
                    if (edges[i].from != edges[i-1].from)
                        currLabelSources++;
                }
                else {
                    if (currLabelCount > 2) {
#ifdef DEBUG
                        fprintf(stderr,"summary: label: %4d  count: %4d  sources: %4d  FILE: %4d\n",
                                edges[i-1].from_label, currLabelCount, currLabelSources, edges[i-1].from);
#else
                        if (currLabelSources==1)
                            fprintf(stdout,"T NODE: %4d\n", edges[i-1].from);
#endif
                    }
                    idx = i;
                    currLabelCount = 1;
                    currLabelSources = 1;
                }
            }
            if ((idx <edgecount) && (currLabelCount > 2)) {
#ifdef DEBUG
                fprintf(stderr,"summary: label: %4d  count: %4d  sources: %4d  FILE: %4d\n",
                        edges[idx].from_label, currLabelCount, currLabelSources, edges[idx].from);
#else
                if (currLabelSources==1)
                    fprintf(stdout,"T NODE: %4d\n", edges[idx].from);
#endif
            }      


            free(edges);

        }
        else {
            fprintf(stdout,"No edges\n");
        }

        free(marked);
    }
#endif

    fclose(fp);

    /* Step 5: Clean up */
    free(infilename);
    free(outfilename);
	free(seedsfile);
    free(graph_type);
    free(alg_type);
    free_graph(g);
    free(g);
    free(membership);
	free(seeds);

    return 0;
}

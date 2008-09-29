#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "utils.h"
#include "sprng.h"

int main(int argc, char** argv) {

    dyn_graph_t* G;
    long i, j;
    long n, m;
    attr_id_t offset, avg_degree;
    double *params, a, b, c, d;
    double av, bv, cv, dv, S, p;
    long num_batches, num_edges_per_batch;
    double ins_del_ratio;
    int SCALE;
    double var;
    long step;
    char* filename;
    long permute_vertices;
    attr_id_t* permV, tmpVal;
    attr_id_t* tmpArr;
    long u, v;
    int* stream, seed;
    attr_id_t *src;
    attr_id_t *dest;
    attr_id_t *degree;
    attr_id_t *degree_hist;
    attr_id_t *mem_chunk;
    double el_time;

    /* For detecting duplicates */
    dyn_array_t* adj;
    long num_self_edges = 0, num_dups = 0, num_mallocs = 0, num_init_mallocs = 0;
        
    if (argc != 2) {
        fprintf(stderr, "Usage: ./test_dyn_ds <config_file_name>\n");
        exit(-1);
    }
  
    filename = (char *) malloc(1000*sizeof(char)); 
    strcpy(filename, argv[1]);

    G = (dyn_graph_t *) calloc(1, sizeof(dyn_graph_t));

    params = (double *) malloc(10*sizeof(double));
    read_dyn_test_config_file(G, filename, params);
    a = params[0];
    b = params[1];
    c = params[2];
    assert(a+b+c < 1);
    d = 1  - (a+b+c);
    permute_vertices = (long) params[3];
    num_batches = (long) params[4]; 
    ins_del_ratio = params[5];

    n = G->n;
    m = G->m;
    num_edges_per_batch = m/num_batches;

    num_self_edges = 0;
    num_dups = 0;

    mem_chunk = (attr_id_t *) malloc(2 * m  * sizeof(attr_id_t));
    adj = (dyn_array_t *) calloc(n, sizeof(dyn_array_t));
     
    src = (attr_id_t *) malloc (m * sizeof(attr_id_t));
    dest = (attr_id_t *) malloc(m * sizeof(attr_id_t));
    degree = (attr_id_t *) calloc(n, sizeof(attr_id_t));
    degree_hist = (attr_id_t *) calloc((n/10), sizeof(attr_id_t));

    assert(src != NULL);
    assert(dest != NULL);
    assert(degree != NULL);
    assert(degree_hist != NULL);
    assert(mem_chunk != NULL);
    assert(adj != NULL);

    offset = 0;
    avg_degree = (2*m)/n;
    for (i=0; i<n; i++) {
        adj[i].vals = &mem_chunk[offset];
        adj[i].max_size = avg_degree;
        offset += avg_degree;
    }

    /* Initialize RNG stream */
    seed = 2387;
    stream = init_sprng(0, 0, 1, seed, SPRNG_DEFAULT);
    
    SCALE = log2(n);
    fprintf(stderr, "Scale: %d\n", SCALE);
    fprintf(stderr, "Avg degree: %ld\n", avg_degree);

    var = 0.1;
    
    /* Generate edges */
    for (i=0; i<m; i++) {

        u = 1;
        v = 1;
        step = n/2;

        av = a;
        bv = b;
        cv = c;
        dv = d;

        p = sprng(stream);
        if (p < av) {
            /* Do nothing */
        } else if ((p >= av) && (p < av+bv)) {
            v += step;
        } else if ((p >= av+bv) && (p < av+bv+cv)) {
            u += step;
        } else {
            u += step;
            v += step;
        }
        
        for (j=1; j<SCALE; j++) {
            step = step/2;

            /* Vary a,b,c,d by up to 10% */
            av *= 0.95 + var * sprng(stream);
            bv *= 0.95 + var * sprng(stream);
            cv *= 0.95 + var * sprng(stream);
            dv *= 0.95 + var * sprng(stream);

            S = av + bv + cv + dv;
            av = av/S;
            bv = bv/S;
            cv = cv/S;
            dv = dv/S;
            
            /* Choose partition */
            p = sprng(stream);
            if (p < av) {
                /* Do nothing */
            } else if ((p >= av) && (p < av+bv)) {
                v += step;
            } else if ((p >= av+bv) && (p < av+bv+cv)) {
                u += step;
            } else {
                u += step;
                v += step;
            }
        }
       
        if (u == v) {
            num_self_edges++;
            i--;
            continue;
        }

        u--;
        v--;
        if (adj[u].count < adj[v].count) {
            /* Check for duplicate in u's adjacencies */
            for (j = 0; j < adj[u].count; j++) {
                if (adj[u].vals[j] == v) {
                    i--;
                    num_dups++;
                    break;
                }
            }
            if (j < adj[u].count)
                continue;
        } else {
            /* Check for duplicate in v's adjacencies */
            for (j = 0; j < adj[v].count; j++) {
                if (adj[v].vals[j] == u) {
                    i--;
                    num_dups++;
                    break;
                }
            }
            if (j < adj[v].count)
                continue;
        }
        
        /* fprintf(stderr, "%ld %ld ", u, v); */
        if (adj[u].count == adj[u].max_size) {
            num_mallocs++;
            if (adj[u].max_size == avg_degree) {
                num_init_mallocs++;
                tmpArr = (attr_id_t *) malloc(4 * adj[u].max_size *
                     sizeof(attr_id_t));
                assert(tmpArr != NULL);
                memcpy(tmpArr, adj[u].vals, adj[u].max_size *
                    sizeof(attr_id_t));
                adj[u].vals = NULL;
                adj[u].vals = tmpArr;
                adj[u].max_size = adj[u].max_size * 4;
            } else {
                tmpArr = (attr_id_t *) malloc(4 * adj[u].max_size *
                     sizeof(attr_id_t));
                assert(tmpArr != NULL);
                memcpy(tmpArr, adj[u].vals, adj[u].max_size *
                    sizeof(attr_id_t));
                free(adj[u].vals);
                adj[u].max_size = adj[u].max_size * 4;
                adj[u].vals = tmpArr;
            }
        }
        /* fprintf(stderr, "[u %ld %ld] ", adj[u].count, adj[u].max_size); */
        adj[u].vals[adj[u].count++] = v;    
        if (adj[v].count == adj[v].max_size) {
            num_mallocs++;
            if (adj[v].max_size == avg_degree) {
                num_init_mallocs++;
                tmpArr = (attr_id_t *) malloc(4 * adj[v].max_size *
                     sizeof(attr_id_t));
                assert(tmpArr != NULL);
                memcpy(tmpArr, adj[v].vals, adj[v].max_size *
                    sizeof(attr_id_t));
                adj[v].vals = NULL;
                adj[v].vals = tmpArr;
                adj[v].max_size = adj[v].max_size * 4;
            } else {
                tmpArr = (attr_id_t *) malloc(4 * adj[v].max_size *
                     sizeof(attr_id_t));
                assert(tmpArr != NULL);
                memcpy(tmpArr, adj[v].vals, adj[v].max_size *
                    sizeof(attr_id_t));
                free(adj[v].vals);
                adj[v].max_size = adj[v].max_size * 4;
                adj[v].vals = tmpArr;
            }
        }
        /* fprintf(stderr, "[v %ld %ld] ", adj[v].count, adj[v].max_size); */
        adj[v].vals[adj[v].count++] = u;    
        
        src[i] = u;
        dest[i] = v;
    }

    fprintf(stderr, "num self edges: %ld, num dups: %ld\n", num_self_edges,
            num_dups);
    fprintf(stderr, "num mallocs: %ld, num init mallocs: %ld\n", num_mallocs,
            num_init_mallocs);
    
    for (j = 0; j<n; j++) {
        if (adj[j].max_size != avg_degree)
            free(adj[j].vals);
    }

    free(mem_chunk);
    free(adj);

    if (permute_vertices) {

        permV = (attr_id_t *) malloc(n*sizeof(attr_id_t));
        assert(permV != NULL);

        for (i=0; i<n; i++) {
            permV[i] = i;
        }

        for (i=0; i<n; i++) {
            j = n * sprng(stream);
            tmpVal = permV[i];
            permV[i] = permV[j];
            permV[j] = tmpVal;
        }

        for (i=0; i<m; i++) {
            src[i]  = permV[src[i]];
            dest[i] = permV[dest[i]];
       }

    }

    for (i=0; i<m; i++) {
        degree[src[i]]++;
        /* fprintf(stderr, "[%d %d] ", src[i], dest[i]); */
    }

    for (i=0; i<n; i++) {
        degree_hist[degree[i]]++;
    }

    for (i=0; i<n/10; i++) {
        if (degree_hist[i] != 0) {
            fprintf(stderr, "%ld %ld\n", i, degree_hist[i]);
        }
    }

    /* Case 1: Oracle, we know bucket sizes */
    for (i=1; i<n; i++) {
        degree[i] += degree[i-1];
    }
    mem_chunk = (attr_id_t *) malloc(2 * m  * sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(n, sizeof(dyn_array_t));

    G->adj[0].vals = mem_chunk;
    G->adj[0].max_size = degree[0];

    for (i=1; i<G->n; i++) {
        G->adj[i].vals = mem_chunk + degree[i];
        G->adj[i].max_size = degree[i]-degree[i-1];
    }
    
#ifdef _OPENMP
#pragma omp parallel
{
#endif
    long iter;
    attr_id_t uid, vid, edge_id;
    
#ifdef _OPENMP
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
#pragma omp single
   fprintf(stderr, "No. of threads: %d\n", nthreads);
#pragma omp barrier
    double elt = omp_get_wtime();
#else
    int tid = 0;
    int nthreads = 1;
    double elt = get_seconds();
#endif

#ifdef _OPENMP
#pragma omp for nowait
#endif
    for (iter = 0; iter < m; iter++) {
        uid = src[iter];
        vid = dest[iter];

#ifdef _OPENMP        
#pragma omp atomic
#endif
        edge_id = G->adj[uid].count++;
        G->adj[uid].vals[edge_id] = vid;
    }

#ifdef _OPENMP
    elt = omp_get_wtime() - elt;
    fprintf(stderr, "Thread: %tid, time: %lf sec\n", elt); 
#else
    elt = get_seconds() - elt;
    fprintf(stderr, "Time taken: %lf sec\n", elt);
#endif

#ifdef _OPENMP
}
#endif
    
    /* Case 2: We resize buckets inside code, locking */

    /* Case 3: Oracle, no locking */ 

    /* Case 4: Resize, no locking */

    /* Case 5: Oracle, bunched high-degree vertices */

    /* Case 6: Spanning Forest construction */

    /* Case 7: Deletions, dyn array */

    /* Case 8: Deletions, sorted dyn array */

    /* Case 9: Insertions, Oracle, binary heap for high degree vertices */

    /* Case 10: Deletions, Oracle, binary heap for high degree vertices */


    free(filename);
    free(params);
    free(src);
    free(dest);
    free(degree);
    free(degree_hist);
    free(permV);
    
    
    free(G);
    return 0;
}


void read_dyn_test_config_file(dyn_graph_t* G, char* configfile, double* params) {

    /* read parameters from config file */
    FILE *fp;
    char line[128], var[32];
    double val;

    fp = fopen(configfile,"r");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open config file: %s\n",configfile);
        exit(-1);
    }

    while (fgets(line, sizeof (line), fp) != NULL) {
        sscanf(line, "%s %lf", var, &val);
        if (*var == '#') continue;  /* comment */
        if (strcmp(var, "n") == 0) {
            G->n = (long) val;
            assert(G->n > 0);
        } else if (strcmp(var, "m") == 0) {
            G->m = (long) val;
            assert(G->m > 0);
        } else if (strcmp(var, "a") == 0) {
            params[0] = val;
            assert((val > 0) && (val < 1));
        } else if (strcmp(var, "b") == 0) {
            params[1] = val;
            assert((val > 0) && (val < 1));
        } else if (strcmp(var, "c") == 0) {
            params[2] = val;
            assert((val > 0) && (val < 1)); 
        } else if (strcmp(var, "permute_vertices") == 0) {
            params[3] = val;
        } else if (strcmp(var, "num_batches") == 0) {
            params[4] = val;
        } else if (strcmp(var, "ins_del_ratio") == 0) {
            params[5] = val;
        } else {
            fprintf(stderr,"Unknown parameter: %s\n", line);
        }
    }
    fclose(fp);

}


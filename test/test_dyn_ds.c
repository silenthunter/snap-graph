#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "utils.h"
#include "sprng.h"

int main(int argc, char** argv) {

    dyn_graph_t* G;
    double *params;
    char* filename;
    attr_id_t *src, *dest, *degree;

    if (argc != 2) {
        fprintf(stderr, "Usage: ./test_dyn_ds <config_file_name>\n");
        exit(-1);
    }

    filename = (char *) malloc(1000*sizeof(char)); 
    strcpy(filename, argv[1]);

    G = (dyn_graph_t *) calloc(1, sizeof(dyn_graph_t));
    params = (double *) malloc(10*sizeof(double));
    
    read_dyn_test_config_file(G, filename, params);
    
    src = (attr_id_t *) malloc (G->m * sizeof(attr_id_t));
    dest = (attr_id_t *) malloc(G->m * sizeof(attr_id_t));
    degree = (attr_id_t *) calloc(G->n, sizeof(attr_id_t));

    assert(src != NULL);
    assert(dest != NULL);
    assert(degree != NULL);

    /* Generate edges */
    par_gen_RMAT_edges(G, params, src, dest, degree);

    /* Case 1: Oracle, we know bucket sizes */
    dyn_ds_init_nomalloc(G, src, dest, degree);

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
    free(G);

    
    return 0;
}

void par_gen_RMAT_edges(dyn_graph_t* G, double* params, attr_id_t* src,
        attr_id_t* dest, attr_id_t* degree) {

    attr_id_t *permV;
    dyn_array_t* adj;
    attr_id_t *degree_hist;
    attr_id_t* edgeCount;
    attr_id_t *mem_chunk;
    long num_self_edges, num_dups;
    long num_mallocs, num_init_mallocs;
    omp_lock_t* vLock;
    attr_id_t edges_added, curr_num_edges_added;
    num_dups = num_self_edges = num_mallocs = num_init_mallocs = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:num_dups) reduction(+:num_self_edges) \
   reduction(+:num_mallocs) reduction(+:num_init_mallocs) num_threads(24)
#endif
{
    double el_time;
    long n, m;
    attr_id_t offset, avg_degree;
    long num_batches, num_edges_per_batch;
    
    long i, j;
    double a, b, c, d;
    double av, bv, cv, dv, S, p;
    double ins_del_ratio;
    int SCALE;
    double var;
    long step;
    long permute_vertices;
    attr_id_t tmpVal;
    attr_id_t* tmpArr;
    long u, v;
    int* stream, seed;
    attr_id_t brk;
    int tid, nthreads;
    attr_id_t *psrc, *pdest;
    attr_id_t pEdgeCount;
    attr_id_t num_edges_to_add;

#ifdef _OPENMP
    nthreads = omp_get_num_threads();
    tid = omp_get_thread_num(); 
#else
    nthreads = 1;
    tid = 0;
#endif

    if (tid == 0) 
        el_time = get_seconds();
   
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
    // num_edges_per_batch = m/num_batches;

    num_self_edges = 0;
    num_dups = 0;

    if (tid == 0) {
        mem_chunk = (attr_id_t *) malloc(2 * m  * sizeof(attr_id_t));
        adj = (dyn_array_t *) calloc(n, sizeof(dyn_array_t));
        degree_hist = (attr_id_t *) calloc((n/10), sizeof(attr_id_t));
        vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
        edgeCount = (attr_id_t *) malloc((nthreads+1)*sizeof(attr_id_t));
        
        assert(degree_hist != NULL);
        assert(mem_chunk != NULL);
        assert(adj != NULL);
        assert(vLock != NULL);
        assert(edgeCount != NULL);
    }
    
    avg_degree = (2*m)/n;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
    for (i=0; i<n; i++) {
        omp_init_lock(&vLock[i]);
        adj[i].vals = &mem_chunk[i*avg_degree];
        adj[i].max_size = avg_degree;
    }

    if (tid == 0) {
        el_time = get_seconds() - el_time;
        fprintf(stderr, "Init time: %lf\n", el_time);
        el_time = get_seconds();
        fprintf(stderr, "Generating edges ... ");
        edges_added = 0; 
    }
 
    /* Initialize RNG stream */
    seed = 2387;
    stream = init_sprng(0, tid, nthreads, seed, SPRNG_DEFAULT);
    SCALE = log2(n);
    var = 0.1;

    num_edges_to_add = 1<<18;

    psrc = (attr_id_t *) malloc(num_edges_to_add * sizeof(attr_id_t));
    pdest = (attr_id_t *) malloc(num_edges_to_add * sizeof(attr_id_t));

#ifdef _OPENMP
#pragma omp barrier
#endif

    /* Generate edges */
    while (edges_added < m) {

        if (num_edges_to_add > m - edges_added) {
            num_edges_to_add = m - edges_added;
        }

        pEdgeCount = 0;
#ifdef _OPENMP
#pragma omp barrier
#pragma omp for nowait
#endif
        for (i=0; i<num_edges_to_add; i++) {

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
                continue;
            }

            u--;
            v--;
            brk = 0;
            if (adj[u].count < adj[v].count) {
                omp_set_lock(&vLock[u]);
                /* Check for duplicate in u's adjacencies */
                for (j = 0; j < adj[u].count; j++) {
                    if (adj[u].vals[j] == v) {
                        num_dups++;
                        break;
                    }
                }
                if (j < adj[u].count)
                    brk = 1;
                omp_unset_lock(&vLock[u]);
                if (brk)
                    continue;
           
            } else {
                omp_set_lock(&vLock[v]);
                /* Check for duplicate in v's adjacencies */
                for (j = 0; j < adj[v].count; j++) {
                    if (adj[v].vals[j] == u) {
                        num_dups++;
                        break;
                    }
                }
                if (j < adj[v].count)
                    brk = 1;
                omp_unset_lock(&vLock[v]);
                if (brk)
                    continue;
            }
        
            /* fprintf(stderr, "%ld %ld ", u, v); */
            omp_set_lock(&vLock[u]);
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
            omp_unset_lock(&vLock[u]);

            omp_set_lock(&vLock[v]);
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
            omp_unset_lock(&vLock[v]);
            psrc[pEdgeCount] = u;
            pdest[pEdgeCount] = v;
            pEdgeCount++;
            // if (i % 10000 == 0) {
            //     fprintf(stderr, "%ld ", i);
            // }
        }

        edgeCount[tid+1] = pEdgeCount;
#ifdef _OPENMP
#pragma omp barrier
#endif
        if (tid == 0) {
        
            edgeCount[0] = 0; 
            for (j=1; j<nthreads; j++) {
                edgeCount[j+1] += edgeCount[j];
            }
        }

#ifdef _OPENMP
#pragma omp barrier
#endif
        /* Write partial results to src and dest */
        int offset = edges_added + edgeCount[tid];
        for (j = 0; j<pEdgeCount; j++) {
            src[offset+j] = psrc[j];
            dest[offset+j] = pdest[j];
        }

#ifdef _OPENMP
#pragma omp barrier
#endif
        if (tid == 0) {
            edges_added += edgeCount[nthreads];
            fprintf(stderr, "%ld ", edges_added);
        }
#ifdef _OPENMP
#pragma omp barrier
#endif
    }

    free(psrc);
    free(pdest);

    if (tid == 0) {
        fprintf(stderr, "\n");
        el_time = get_seconds() - el_time;
        fprintf(stderr, "Edge generation time: %lf sec\n", el_time);
    }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
    for (j = 0; j<n; j++) {
        if (adj[j].max_size != avg_degree)
            free(adj[j].vals);
    }

#ifdef _OPENMP
#pragma omp barrier
#endif
    if (tid == 0) {
        free(mem_chunk);
        free(adj);
        free(edgeCount);
        el_time = get_seconds();
    }

#ifdef _OPENMP
#pragma omp barrier
#endif

    if (permute_vertices) {

        if (tid == 0) {
            permV = (attr_id_t *) malloc(n*sizeof(attr_id_t));
            assert(permV != NULL);
        }
 
#ifdef _OPENMP
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp for
#endif       
        for (i=0; i<n; i++) {
            permV[i] = i;
        }

#ifdef _OPENMP
#pragma omp for
#endif
        for (i=0; i<n; i++) {
            j = n * sprng(stream);
            if (i != j) {
                omp_set_lock(&vLock[i]);
                omp_set_lock(&vLock[j]);
                tmpVal = permV[i];
                permV[i] = permV[j];
                permV[j] = tmpVal;
                omp_unset_lock(&vLock[j]);
                omp_unset_lock(&vLock[i]);
            }
        }

#ifdef _OPENMP
#pragma omp for
#endif
        for (i=0; i<m; i++) {
            // fprintf(stderr, "[%ld %ld] ", src[i], dest[i]);
            src[i]  = permV[src[i]];
            dest[i] = permV[dest[i]];
        }

        if (tid == 0) {
            el_time = get_seconds() - el_time;
            fprintf(stderr, "Permutation time: %lf sec\n", el_time);
            el_time = get_seconds();
        }
    }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
#endif
    for (i=0; i<m; i++) {
#ifdef _OPENMP
#pragma omp atomic
#endif
        degree[src[i]]++;
    }

#ifdef _OPENMP
#pragma omp for
#endif
    for (i=0; i<n; i++) {
#ifdef _OPENMP
#pragma omp atomic
#endif
        degree_hist[degree[i]]++;
    }

    if (tid == 0) {
#if 0
        for (i=0; i<n/10; i++) {
            if (degree_hist[i] != 0) {
                 fprintf(stderr, "[%ld %ld] ", i, degree_hist[i]);
            }
        }
#endif
        el_time = get_seconds() - el_time;
        fprintf(stderr, "\nDegree hist time: %lf sec\n", el_time);
    }

#ifdef _OPENMP
#pragma omp for 
#endif
    for (int i=0; i<n; i++) {
        omp_destroy_lock(&vLock[i]);
    }

}

    fprintf(stderr, "num self edges: %ld, num dups: %ld\n", num_self_edges,
            num_dups);
    fprintf(stderr, "num mallocs: %ld, num init mallocs: %ld\n", num_mallocs,
            num_init_mallocs);

    free(vLock);
    free(degree_hist);
    free(permV);
    
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

void dyn_ds_init_nomalloc(dyn_graph_t* G, attr_id_t* src, attr_id_t* dest,
        attr_id_t* degree) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
#endif
    attr_id_t* mem_chunk;

    el_time = get_seconds();
    for (int i=1; i<G->n; i++) {
        degree[i] += degree[i-1];
    }
    mem_chunk = (attr_id_t *) malloc(G->m  * sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    G->adj[0].vals = mem_chunk;
    G->adj[0].max_size = degree[0];

    for (int i=1; i<G->n; i++) {
        G->adj[i].vals = mem_chunk + degree[i];
        G->adj[i].max_size = degree[i]-degree[i-1];
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(24)
#endif
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
#pragma omp parallel num_threads(1)
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    double elt;
    int tid, nthreads;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
#pragma omp single
    fprintf(stderr, "No. of threads: %d\n", nthreads);
#pragma omp barrier
    elt = omp_get_wtime();
#else
    tid = 0;
    nthreads = 1;
    elt = get_seconds();
#endif

#ifdef _OPENMP
#pragma omp for
#endif
    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
        
#ifdef _OPENMP        
        omp_set_lock(&vLock[u]);
        eid = G->adj[u].count++;
        omp_unset_lock(&vLock[u]);
#else
        eid = G->adj[u].count++;
#endif
        fprintf(stderr, "[%ld %ld %ld] ", u, v, eid);
        G->adj[u].vals[eid] = v;
    }

#ifdef _OPENMP
    elt = omp_get_wtime() - elt;
    fprintf(stderr, "Thread: %ld, time: %lf sec\n", tid, elt); 
#pragma omp barrier
#else
    elt = get_seconds() - elt;
    fprintf(stderr, "Time taken: %lf sec\n", elt);
#endif
}

#ifdef _OPENMP
#pragma omp parallel for
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif
    
    free(mem_chunk);
    free(G->adj);
#ifdef _OPENMP
    free(vLock);
#endif

}

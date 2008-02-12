#include "defs.h"

void SSSP(graph *G, int s, int* weights, int* dist) {

    long n, m;
    long i, j, u, v;
    long tent_d;
    binaryHeap* H;
    
    n = G->n;
    m = G->m;
    
    for (i=0; i<n; i++) {
        dist[i] = INFTY;
    }
    
    dist[s] = 0;
    
    H = (binaryHeap *) malloc(sizeof(binaryHeap));

    initHeap(H, n);
    insertKey(H, s, 0);    
    
    while (H->count != 0) {
        
        u = extractMin(H);
        // fprintf(stderr, "s: %d, d[s]: %d\n", u, dist[u]);
        for (i = G->numEdges[u]; i<G->numEdges[u+1]; i++) {
            v = G->endV[i];
            tent_d = dist[u] + weights[i];

            if (tent_d < dist[v]) {
                if (dist[v] == INFTY) {
                    dist[v] = tent_d;
                    insertKey(H, v, dist[v]);
                    // fprintf(stderr, "%d, weight %d inserted\n", v, dist[v]);
                } else {
                    dist[v] = tent_d;
                    decreaseKey(H, v, dist[v]);
                    // fprintf(stderr, "%d, weight %d decrease_key\n", v, dist[v]);
                }
            }
        }
        /*
        for (i=1; i<=H->count; i++) {
            fprintf(stderr, "<%d %d %d> ", H->nodes[i].key, H->nodes[i].value,
                    H->nodePosition[H->nodes[i].key]);
        }
        fprintf(stderr, "\n");
        */
    }

    freeHeap(H);
}

void BFS_UY(graph *G, int s) {

    /* Dv: set of distinguished vertices */
    int *Dv;
    int Dvcount;
    int *Dvmask;
    int added_count;
    
    /* Dvdist: distance from source vertex to distinguished 
       vertices, determined after constructing the distinguished 
       vertices graph in phase 2 and running SSSP from s*/
    int *Dvdist;
    
    /* Final output: distance from source to all other vertices */
    int *dist;

    /* Each vertex maintains a list of distibguished vertices that 
       it is reachable from, and the distance to the dist. vertex */
    dynArray **reachableD_v;
    dynArray **reachableD_dist;

    /* the graph constructed with distinguished vertices */
    graph *DG;
    long *DGnumEdges;
    int *DGendV;
    int *DGweights;
   
    double t1;

#ifdef _OPENMP    
    omp_lock_t *vLock;
#endif
   
    t1 = get_seconds();

    /* Phase 1 */
    /* create a set of random distinguished vertices */
    
    Dvcount = ceil(sqrt(G->n));
    Dv = (int *) malloc(Dvcount * sizeof(int));
    assert(Dv != NULL);

    Dvmask = (int *) malloc(G->n * sizeof(int));
    assert(Dvmask != NULL);

    for (int i=0; i<G->n; i++) {
       Dvmask[i] = -1;
    }
    fprintf(stderr, "Dvcount: %d\n", Dvcount);

    Dv[0] = s;
    Dvmask[s] = 0;
    added_count = 1; 
    
    srand48(23423);
    while (added_count < Dvcount) {
        long randSrc = lrand48() % G->n;
        if (Dvmask[randSrc] == -1) {
            Dv[added_count] = (int) randSrc;
            Dvmask[randSrc] = added_count;
            added_count++;
        }
    }

    reachableD_v = (dynArray **) malloc(Dvcount * sizeof(dynArray *));
    assert(reachableD_v != NULL); 
    reachableD_dist = (dynArray **) malloc(Dvcount * sizeof(dynArray *));
    assert(reachableD_dist != NULL);

    for (int i=0; i<Dvcount; i++) {
        reachableD_v[i] = (dynArray *) malloc(sizeof(dynArray));
        reachableD_dist[i] = (dynArray *) malloc(sizeof(dynArray));
        assert(reachableD_v[i] != NULL);
        assert(reachableD_dist[i] != NULL);
        initDynArray(reachableD_v[i]);
        initDynArray(reachableD_dist[i]);
    }

    /* Initialize DG */
    DG = (graph *) malloc(sizeof(graph));
    assert(DG != NULL);

    DG->n = Dvcount;
    DGnumEdges = (long *) calloc((Dvcount+1), sizeof(long));
    assert(DGnumEdges != NULL);
    Dvdist = (int *) calloc(Dvcount, sizeof(int));
    assert(Dvdist != NULL);

    dist = (int *) malloc(G->n * sizeof(int));
    assert(dist != NULL);

#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n * sizeof(omp_lock_t));
    assert(vLock != NULL);
#endif

    t1 = get_seconds() - t1;
    fprintf(stderr, "Phase 1 time: %lf\n", t1);

#ifdef _OPENMP
#pragma omp parallel
#endif
{

    int tid;
    long i, j, u, v, w;
    long maxPathLen;
    dynArray **adjD_v;
    dynArray **adjD_d;
    int *S;
    long start, end;
    int *d;
    long n;
#ifdef _OPENMP    
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    n = G->n;
    
    /* Phase 1: path-limited BFSes from the distinguished vertices */
   
    S = (int *) malloc(n*sizeof(int));
    d = (int *) malloc(n*sizeof(int));
    assert(S != NULL);
    assert(d != NULL);
    
    for (i=0; i<n; i++) {
        d[i] = -1;
    }

    adjD_v = (dynArray **) malloc(n * sizeof(dynArray *));
    adjD_d = (dynArray **) malloc(n * sizeof(dynArray *));

    assert(adjD_v != NULL);
    assert(adjD_d != NULL);

    for (i=0; i<n; i++) {
        adjD_v[i] = (dynArray *) malloc(sizeof(dynArray));
        adjD_d[i] = (dynArray *) malloc(sizeof(dynArray));
        assert(adjD_v[i] != NULL);
        assert(adjD_d[i] != NULL);
    
        initDynArray(adjD_v[i]);
        initDynArray(adjD_d[i]);
    }
    
    maxPathLen = 6; // n/Dvcount + 1;
    fprintf(stderr, "maxPathLen: %ld\n", maxPathLen);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (i=0; i<Dvcount; i++) {
        int src = Dv[i];
        // fprintf(stderr, "Source: %d, tid %d\n", src, tid);
        if (G->numEdges[src+1] - G->numEdges[src] == 0)
            continue;

        d[src] = 0;
        start = 0;
        end = 1;
        S[0] = src; 

        while (end - start > 0) {
            
            v = S[start];
            
            if (d[v] >= maxPathLen) {
                break;
            }

            for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {
                w = G->endV[j];

                if (v != w) {

                    if (d[w] < 0) {
                        S[end++] = (int) w;
                        d[w] = d[v] + 1;
                        dynArrayInsert(adjD_v[w], src);
                        dynArrayInsert(adjD_d[w], d[w]);

                        if (Dvmask[w] != -1) {
                            dynArrayInsert(reachableD_v[Dvmask[src]], (int) w);
                            dynArrayInsert(reachableD_dist[Dvmask[src]], d[w]);
                        }
                    }

                }

            }

            start++;
        }

        /* clear all vars */
        for (j=0; j<end; j++) {
            v = S[j];   
            d[v] = -1;    
        }
    }
 
    
    /* Phase 2: construct a graph from the distinguished nodes */
    if (tid == 0) {
        DGnumEdges[0] = 0;
        for (i=1; i<Dvcount; i++) {
            DGnumEdges[i] = DGnumEdges[i-1] + reachableD_v[i-1]->count;
            // fprintf(stderr, "[%ld %ld] ", i-1, reachableD_v[i-1]->count); 
        }
        DG->m = DGnumEdges[Dvcount-1];
        fprintf(stderr, "No. of edges in distinguished graph: %ld\n", DG->m);
        DGendV = (int *) malloc(DG->m * sizeof(int));
        assert(DGendV != NULL);

        DGweights = (int *) malloc(DG->m * sizeof(int));
        assert(DGweights != NULL);
    }
#ifdef _OPENMP
    #pragma omp barrier
#endif

#ifdef _OPENMP    
    #pragma omp for schedule(dynamic)
#endif
    for (i=0; i<Dvcount; i++) {
        for (j=DGnumEdges[i]; j<DGnumEdges[i+1]; j++) {
            DGendV[j] = Dvmask[reachableD_v[i]->vals[j-DGnumEdges[i]]];
            DGweights[j] = reachableD_dist[i]->vals[j-DGnumEdges[i]];
        }
    }

    /* Compute shortest paths */
    if (tid == 0) {
        DG->endV = DGendV;
        DG->numEdges = DGnumEdges;

        for (i=0; i<DG->n; i++) {
            for (j=DG->numEdges[i]; j<DG->numEdges[i+1]; j++) {
                fprintf(stderr, "[%d %d %d] ", i, DG->endV[j], DGweights[j]);
            }
        }
        fprintf(stderr, "done\n"); 
        SSSP(DG, 0, DGweights, Dvdist);
    }

#ifdef _OPENMP
    #pragma omp barrier
#endif

#ifdef _OPENMP    
    #pragma omp for schedule(static)
    for (i=0; i<n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
 
#ifdef _OPENMP    
    #pragma omp for schedule(static)
#endif
    for (i=0; i<n; i++) {    
        dist[i] = INFTY;
    }

    /* Phase 3: determine the distance from source to all other vertices */
    if (tid == 0) {
        dist[s] = 0;
        Dvdist[Dvmask[s]] = 0;
    }
#ifdef _OPENMP  
    #pragma omp barrier
#endif
    
    for (i=0; i<n; i++) {
        if (Dvmask[i] != -1) {
            dist[i] = Dvdist[Dvmask[i]];
            continue;
        }
        // fprintf(stderr, "count: %d, thread %d\n", adjD_v[i]->count, tid);
        for (j=0; j<adjD_v[i]->count; j++) {
            u = adjD_v[i]->vals[j];
            int dt = Dvdist[Dvmask[u]] + adjD_d[i]->vals[j];
            if (dt < dist[i]) {
#ifdef _OPENMP
                omp_set_lock(&vLock[i]);
#endif
                if (dt < dist[i]) {
                    dist[i] = dt;
                }
#ifdef _OPENMP
                omp_unset_lock(&vLock[i]);
#endif
            }
        } 
    }

#ifdef _OPENMP
    #pragma omp barrier
#endif    
    
#ifdef _OPENMP
    #pragma omp for schedule(static)
    for (i=0; i<n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif

    free(S);
    free(d);

    for (i=0; i<n; i++) {
        freeDynArray(adjD_v[i]);
        free(adjD_v[i]);
        freeDynArray(adjD_d[i]);
        free(adjD_d[i]);
    }
   
    free(adjD_v);
    free(adjD_d);
}    

    free(Dv);
    free(Dvmask);
    free(Dvdist);

    for (int i=0; i<G->n; i++) {
        fprintf(stderr, "%d %d\n", i, dist[i]);
    }

    free(dist);

    for (int i=0; i<Dvcount; i++) {
        freeDynArray(reachableD_v[i]);
        free(reachableD_v[i]);
        freeDynArray(reachableD_dist[i]);
        free(reachableD_dist[i]);
    }
    free(reachableD_v);
    free(reachableD_dist);

    free(DGweights);
    free(DG->numEdges);
    free(DG->endV);
    free(DG);
#ifdef _OPENMP
    free(vLock);
#endif
}

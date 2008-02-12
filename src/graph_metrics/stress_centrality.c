#include "graph_defs.h"
#include "graph_metrics.h"
#include "utils.h"

void stress_centrality_simple(graph_t* G, long* SC) {
#if 0
    long i, j, k, n, chunkSize;
    long tid, nthreads;
    omp_lock_t* vLock;
    long* S;          /* stack of vertices in the order of non-decreasing 
                         distance from s. Also used for implicitly representing 
                         the BFS queue */
    dynArray* P;      /* predecessors of a vertex v on shortest paths from s */
    long* sig;        /* No. of shortest paths */
    long* d;          /* Length of shortest path between every pair */
    long* del;        /* dependency of vertices */
    long start, end;
    long v, w;
    double running_time;


    /* The outer loop is parallelized in this case. Each thread does a BFS 
    and the vertex BC values are incremented atomically */   

#pragma omp parallel \
private(tid, nthreads, S, P, sig, del, d, start, end, v, w) \
private(i, j, k, n)
{
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
     
    n = G->n;

    if (tid == 0) {
        chunkSize = n/nthreads;
        vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
        assert(vLock != NULL);
    }
   
    #pragma omp barrier
    
    #pragma omp for schedule(static, chunkSize)
    for (i=0; i<n; i++) {
        omp_init_lock(&vLock[i]);
    }
    
    S   = (long *) malloc(n*sizeof(long));
    P   = (dynArray *) calloc(n, sizeof(dynArray));
    sig = (long *) malloc(n*sizeof(long));
    d   = (long *) malloc(n*sizeof(long));
    del = (long *) calloc(n, sizeof(long));
   
    assert(S != NULL);
    assert(P != NULL);
    assert(sig != NULL);
    assert(d != NULL);
    assert(del != NULL);
    
    for (i=0; i<n; i++) {
        initDynArray(&P[i]);    
        d[i] = -1;    
    }

    #pragma omp barrier
 
    if (tid == 0) {
        running_time = omp_get_wtime();
    }

    #pragma omp for schedule(dynamic)
    for (i=0; i<n; i++) {
#if VERBOSE_OUTPUT
        if ((i % 1000) == 0)
            fprintf(stdout, "%ld ", i);
#endif
        sig[i] = 1;
        d[i] = 0;

        start = 0;
        S[0] = i;
        end = 1;
        while (end - start > 0) {

            v = S[start];

            for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {

                w = G->endV[j];
                if (v != w) {

                    /* w found for the first time? */ 
                    if (d[w] < 0) {
                        S[end++] = w;
                        d[w] = d[v] + 1;
                        sig[w] = 0;
                    }

                    if (d[w] == d[v] + 1) {
                        sig[w] = sig[w] + sig[v];
                        dynArrayInsert(&P[w], v);
                    }

                }
            }

            start++;

        }
        
        for (j=end-1; j>0; j--) {
            w = S[j];
            for (k = 0; k<P[w].count; k++) {
                v = P[w].vals[k];
                del[v] = del[v] + sig[v]*(1+del[w]);
            }
            /* increment BC value atomically */
            omp_set_lock(&vLock[w]);
            SC[w] += del[w];
            omp_unset_lock(&vLock[w]);
        } 

        for (j=end-1; j>=0; j--) {
            w = S[j];
            d[w] = -1;
            del[w] = 0;
            P[w].count = 0;
        }
    }

    if (tid == 0) { 
        running_time = omp_get_wtime() - running_time;
        fprintf(stderr, "done.\nTime taken: %lf seconds.\n", running_time);
    }
    
    #pragma omp for schedule(static, chunkSize)
    for (i=0; i<n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
     
    free(S);
    for (j=0; j<n; j++)
        freeDynArray(&P[j]);
    free(P);
    free(sig);
    free(d);
    free(del);
    if (tid == 0)
        free(vLock);
    #pragma omp barrier
}    
#endif
}

void stress_centrality_parBFS(graph_t* G, long* SC) {
#if 0
    long i, j, k, n, chunkSize;
    long tid, nthreads;
    omp_lock_t* vLock;
    long *S, *myS;      /* stack of vertices in the order of non-decreasing 
                           distance from s. Also used for implicitly representing 
                           the BFS queue */
    dynArray* P;        /* predecessors of a vertex v on shortest paths from s */
    long* sig;          /* No. of shortest paths */
    long* d;            /* Length of shortest path between every pair */
    long* del;          /* dependency of vertices */
    long *start, *end;
    long *psCount, myCount;
    long v, w, vert;
    int myLock;
    long numV, count, numPhases;
    long MAX_NUM_PHASES;
    double running_time;


#pragma omp parallel \
private(tid, nthreads, v, vert, w, numV, myS, myCount, myLock,numPhases) \
private(i, j, k, n, count)
{
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
     
    numV = G->n;
    n = G->n;

    myCount = 0;
    MAX_NUM_PHASES = 100;
    
    
    if (tid == 0) {
        chunkSize = n/nthreads;
        vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
    }

    #pragma omp barrier

    #pragma omp for schedule(static, chunkSize)
    for (i=0; i<n; i++) {
        omp_init_lock(&vLock[i]);
    }

    if (tid == 0) {
        S   = (long *) malloc(n*sizeof(long));
        P   = (dynArray *) calloc(n, sizeof(dynArray));
        sig = (long *) malloc(n*sizeof(long));
        d   = (long *) malloc(n*sizeof(long));
        del = (long *) calloc(n, sizeof(long));
        start = (long *) malloc(MAX_NUM_PHASES*sizeof(long));
        end = (long *) malloc(MAX_NUM_PHASES*sizeof(long));
        psCount = (long *) malloc((nthreads+1)*sizeof(long));
        
        assert(S != NULL);
        assert(P != NULL);
        assert(sig != NULL);
        assert(d != NULL);
        assert(del != NULL);
        assert(start != NULL);
        assert(end != NULL);
        assert(psCount != NULL);
    }
   
    myS = (long *) malloc(n*sizeof(long));
    myCount = 0;
    
    #pragma omp barrier
     
    #pragma omp for schedule(static, chunkSize) 
    for (i=0; i<n; i++) {
        initDynArray(&P[i]);
        d[i] = -1;
    }

    if (tid == 0) {
        running_time = omp_get_wtime();
    }
    
    #pragma omp barrier
    
    for (i=0; i<numV; i++) {
        if (tid == 0) {
            sig[i] = 1;
            d[i] = 0;
            S[0] = i;
            start[0] = 0;
            end[0] = 1;
        }
        
        count = 1;
        numPhases = 0;
       
        #pragma omp barrier

        while (end[numPhases] - start[numPhases] > 0) {

            myCount = 0;
            #pragma omp barrier
            
            #pragma omp for schedule(dynamic)
            for (vert = start[numPhases]; vert < end[numPhases]; vert++) {
                
                v = S[vert];
                for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {

                    w = G->endV[j];
                    if (v != w) {
                         
                        myLock = omp_test_lock(&vLock[w]);
                           
                        if (myLock) { 
                            /* w found for the first time? */ 
                            if (d[w] == -1) {
                                myS[myCount++] = w;
                                d[w] = d[v] + 1;
                                sig[w] = 0;
                                sig[w] += sig[v];
                                dynArrayInsert(&P[w], v);
                            } else if (d[w] == d[v] + 1) {
                                sig[w] += sig[v];
                                dynArrayInsert(&P[w], v);
                            }
                               
                            omp_unset_lock(&vLock[w]);
                        } else {
                            if ((d[w] == -1) || (d[w] == d[v]+ 1)) {
                                omp_set_lock(&vLock[w]);
                                sig[w] += sig[v];
                                dynArrayInsert(&P[w], v);
                                omp_unset_lock(&vLock[w]);
                            }
                        }
                           
                    }
                }
            }
            
            /* Merge all local stacks for next iteration */    
            numPhases++; 
            
            psCount[tid+1] = myCount;

            #pragma omp barrier
            
            if (tid == 0) {
                start[numPhases] = end[numPhases-1];
                psCount[0] = start[numPhases];
                for(k=1; k<=nthreads; k++) {
                    psCount[k] = psCount[k-1] + psCount[k];
                }
                end[numPhases] = psCount[nthreads];
            }
            #pragma omp barrier
            
            for (k = psCount[tid]; k < psCount[tid+1]; k++) {
                S[k] = myS[k-psCount[tid]];
            } 
            
            #pragma omp barrier
            count = end[numPhases];
        }
     
        numPhases--;
        
        #pragma omp barrier

        while (numPhases > 0) {
            #pragma omp for schedule(dynamic)
            for (j=start[numPhases]; j<end[numPhases]; j++) {
                w = S[j];
                for (k = 0; k<P[w].count; k++) {
                    v = P[w].vals[k];
                    omp_set_lock(&vLock[v]);
                    del[v] = del[v] + sig[v]*(1+del[w]);
                    omp_unset_lock(&vLock[v]);
                }
                SC[w] += del[w];
            }
       
            #pragma omp barrier
            
            numPhases--;

            #pragma omp barrier
            
        }

        if (tid == 0) {
            chunkSize = count/nthreads;
        }
        
        #pragma omp barrier
       
        #pragma omp for schedule(static, chunkSize)
        for (j=0; j<count; j++) {
            w = S[j];
            d[w] = -1;
            del[w] = 0;
            P[w].count = 0;
        }

        #pragma omp barrier

    }

    if (tid == 0) { 
        running_time = omp_get_wtime() - running_time;
        fprintf(stderr, "done.\nTime taken: %lf seconds.\n", running_time);
    }
   

    #pragma omp for schedule(static, chunkSize)
    for (i=0; i<n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
   
    free(myS);
    
    if (tid == 0) { 
        free(S);
        for (j=0; j<n; j++)
            freeDynArray(&P[j]);
        free(P);
        free(sig);
        free(d);
        free(del);
        free(vLock);
        free(start);
        free(end);
        free(psCount);
    }

}    
#endif

}

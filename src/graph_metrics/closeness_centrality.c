#include "graph_defs.h"
#include "graph_metrics.h"
#include "utils.h"

void closeness_centrality_simple(graph_t* G, double* CC) {
#if 0
    long i, j, n;
    long tid;
    long* S;          /* The BFS queue */ 
    long start, end;
    long *d;
    long dist;
    long v, w;
    double running_time;
    
    /* The outer loop is parallelized in this case. Each thread does a BFS 
    from a vertex, and the closeness centrality value is computed */   

#pragma omp parallel \
private(tid, S, start, end, v, w, d, dist) \
private(i, j, n, running_time)
{
    tid = omp_get_thread_num();
     
    n = G->n;

    S   = (long *) malloc(n*sizeof(long));
    assert(S != NULL);
    d   = (long *) malloc(n*sizeof(long));
    assert(d != NULL);
   
    for (i=0; i<n; i++) {
        d[i] = -1;    
    }
    
    if (tid == 0) {
        running_time = omp_get_wtime();
    }

    #pragma omp barrier
    
    #pragma omp for schedule(dynamic)
    for (i=0; i<n; i++) {
#if VERBOSE_OUTPUT
        if ((i % 1000) == 0)
            fprintf(stdout, "%ld ", i);
#endif
        d[i] = 0;
        start = 0;
        S[0] = i;
        end = 1;
        dist = 0;

        while (end - start > 0) {
            v = S[start];
            for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {
                w = G->endV[j];
                if (v != w) {
                    /* w found for the first time? */ 
                    if (d[w] < 0) {
                        S[end++] = w;
                        d[w] = d[v] + 1;
                        dist += d[w]; 
                    }
                }
            }
            start++;
        }
        
        if (dist == 0)
            CC[i] = 0;
        else 
            CC[i] = 1/(double) dist;
        
        for (j=0; j<end; j++) {
            w = S[j];
            d[w] = -1;
        }
    }
    
    #pragma omp barrier

    if (tid == 0) { 
        running_time = omp_get_wtime() - running_time;
        fprintf(stderr, "done.\nTime taken: %lf seconds.\n", running_time);
    }
    
    free(S);
    free(d);
}    
#endif
}

void closeness_centrality_parBFS(graph_t* G, double* CC) {
#if 0
    long i, j, k, n;
    long tid, nthreads;
    omp_lock_t* vLock;
    long *S, *myS;      /* stack of vertices in the order of non-decreasing 
                           distance from s. Also used for implicitly representing 
                           the BFS queue */
    long* d;            /* Length of the shortest path from the source to a vertex */
    long *start, *end;
    long *psCount, myCount;
    long v, w, vert;
    long *dist, dist_sum;
    long chunkSize;
    int myLock;
    long count, numPhases;
    long MAX_NUM_PHASES;
    double running_time;


#pragma omp parallel \
private(tid, nthreads, v, vert, w, myS, myCount, myLock, numPhases) \
private(i, j, k, n, count)
{
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
     
    n = G->n;

    myCount = 0;
    MAX_NUM_PHASES = 100;
    
    if (tid == 0) {
        chunkSize = n/nthreads;
        vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
        dist = (long *) calloc(nthreads, sizeof(long));
        dist_sum = 0;
    }

    #pragma omp barrier

    #pragma omp for schedule(static, chunkSize)
    for (i=0; i<n; i++) {
        omp_init_lock(&vLock[i]);
    }

    if (tid == 0) {
        S   = (long *) malloc(n*sizeof(long));
        d   = (long *) malloc(n*sizeof(long));
        start = (long *) malloc(MAX_NUM_PHASES*sizeof(long));
        end = (long *) malloc(MAX_NUM_PHASES*sizeof(long));
        psCount = (long *) malloc((nthreads+1)*sizeof(long));
        
        assert(S != NULL);
        assert(d != NULL);
        assert(start != NULL);
        assert(end != NULL);
        assert(psCount != NULL);
    }
   
    myS = (long *) malloc(n*sizeof(long));
    assert(myS != 0);

    myCount = 0;
    
    #pragma omp barrier
     
    #pragma omp for schedule(static, chunkSize) 
    for (i=0; i<n; i++) {
        d[i] = -1;
    }

    if (tid == 0) { 
        running_time = omp_get_wtime();
    }

    #pragma omp barrier
    
    for (i=0; i<n; i++) {
        if (tid == 0) {
            d[i] = 0;
            dist_sum = 0;
            S[0] = i;
            start[0] = 0;
            end[0] = 1;
        }
        
        count = 1;
        numPhases = 0;
       
        dist[tid] = 0;

        #pragma omp barrier
        
        while (end[numPhases] - start[numPhases] > 0) {
            myCount = 0;
            
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
                                dist[tid] += d[w];
                            }
                               
                            omp_unset_lock(&vLock[w]);
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

        if (tid == 0) {
            for (j=0; j<nthreads; j++) {
                dist_sum += dist[j];
            }

            CC[i] = 1/(double) dist_sum;
        }
        
        #pragma omp for schedule(static)
        for (j=0; j<count; j++) {
            w = S[j];
            d[w] = -1; 
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
        free(d);
        free(vLock);
        free(start);
        free(end);
        free(psCount);
        free(dist);
    }

}    
#endif
}

#include "graph_defs.h"
#include "graph_kernels.h"

void BFS_sequential(graph_t* G, long src, int* d) {

}

void BFS_parallel_frontier_expansion(graph_t* G, long src, long diameter) {

    attr_id_t* S;
    long *start;
    char* visited;
    long *pSCount;

#ifdef _OPENMP
    omp_lock_t* vLock;
#endif

    long phase_num, numPhases;
    long count;

#ifdef DIAGNOSTIC    
    double elapsed_time = get_seconds();
#endif

    numPhases = diameter + 1;
        
#ifdef _OPENMP 

#pragma omp parallel
{
#endif

    attr_id_t *pS, *pSt;
    long pCount, pS_size;
    long v, w;
    int tid, nthreads;
    
    long j, k, vert, n;

#ifdef _OPENMP    
    long i;
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
#else
    tid = 0;
    nthreads = 1;
#endif
    
    n = G->n;
    
    pS_size = n/nthreads + 1;
    pS = (attr_id_t *) malloc(pS_size*sizeof(attr_id_t));
    assert(pS != NULL);
    
    if (tid == 0) {  
        S = (attr_id_t *) malloc(n*sizeof(attr_id_t));
        visited = (char *) calloc(n, sizeof(char));
        start = (long *) calloc((numPhases+2), sizeof(long));
        pSCount = (long *) malloc((nthreads+1)*sizeof(long));
#ifdef _OPENMP
        vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
#endif
    }

#ifdef _OPENMP    
#pragma omp barrier

#pragma omp for
    for (i=0; i<n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    
#ifdef _OPENMP
#pragma omp barrier
#endif
        
    if (tid == 0) {
        S[0] = src;
        visited[src] = (char) 1;
        count = 1;
        phase_num = 0;
        start[0] = 0;
        start[1] = 1;
    }

        
#ifdef _OPENMP
#pragma omp barrier
#endif

        while (start[phase_num+1] - start[phase_num] > 0) {
        
            pCount = 0;

#ifdef _OPENMP
#pragma omp for
#endif
            for (vert=start[phase_num]; vert<start[phase_num+1]; vert++) {
            
                v = S[vert];
                for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {
                    w = G->endV[j]; 
                    if (v == w)
                        continue;
#ifdef _OPENMP
                    int myLock = omp_test_lock(&vLock[w]);
                    if (myLock) {
#endif
                        if (visited[w] != (char) 1) { 
                            visited[w] = (char) 1;
                            if (pCount == pS_size) {
                                /* Resize pS */
                                pSt = (attr_id_t *)
                                    malloc(2*pS_size*sizeof(attr_id_t));
                                memcpy(pSt, pS, pS_size*sizeof(attr_id_t));
                                free(pS);
                                pS = pSt;
                                pS_size = 2*pS_size;
                            }
                            pS[pCount++] = w;
                        }
#ifdef _OPENMP
                        omp_unset_lock(&vLock[w]);
                    }
#endif
                }
            }
            

#ifdef _OPENMP
#pragma omp barrier
#endif            
            pSCount[tid+1] = pCount;
 
#ifdef _OPENMP
#pragma omp barrier
#endif            
           
            if (tid == 0) {
                pSCount[0] = start[phase_num+1];
                for(k=1; k<=nthreads; k++) {
                    pSCount[k] = pSCount[k-1] + pSCount[k];
                }
                start[phase_num+2] = pSCount[nthreads];
                count = pSCount[nthreads];
                phase_num++;
            }
            
#ifdef _OPENMP
#pragma omp barrier
#endif
            for (k = pSCount[tid]; k < pSCount[tid+1]; k++) {
                S[k] = pS[k-pSCount[tid]];
            } 
            

#ifdef _OPENMP
#pragma omp barrier
#endif
        } /* End of search */

#ifdef DIAGNOSTIC
        if (tid == 0) {
            fprintf(stderr, "Search from vertex %ld, No. of vertices visited: %ld\n", src, count);
        }
#endif

    free(pS);
#ifdef _OPENMP    
#pragma omp barrier

#pragma omp for
    for (i=0; i<n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#pragma omp barrier
#endif

    if (tid == 0) {
        free(S);
        free(start);
        free(visited);
        free(pSCount);
#ifdef _OPENMP
        free(vLock);
#endif
        
    }

#ifdef _OPENMP    
}
#endif

#ifdef DIAGNOSTIC    
    elapsed_time = get_seconds() - elapsed_time;
    fprintf(stderr, "Time taken for BFS: %lf seconds\n", elpased_time);
#endif

}



void BFS_path_limited_search(graph_t* G, long src, long graph_diameter) {

}



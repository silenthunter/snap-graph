#ifndef _GRAPH_DEFS_H
#define _GRAPH_DEFS_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ENABLE_64BIT
    typedef long attr_id_t;
#else
    typedef int attr_id_t;
#endif

#define ARRAY_INIT_SIZE 8

typedef struct {
    attr_id_t* vals;
    long count;
	long max_size;
} dyn_array;

void dyn_array_insert(dyn_array* A, attr_id_t val);
void dyn_array_delete(dyn_array* A, attr_id_t val);
void dyn_array_init(dyn_array* A);
void dyn_array_clear(dyn_array* A);
void dyn_array_free(dyn_array* A);

typedef struct {
    attr_id_t* vals;
    long count;
	long max_size;
} sorted_dyn_array;

void sorted_dyn_array_insert(sorted_dyn_array* A, attr_id_t val);
void sorted_dyn_array_delete(sorted_dyn_array* A, attr_id_t val);
void sorted_dyn_array_init(sorted_dyn_array* A);
void sorted_dyn_array_clear(sorted_dyn_array* A);
void sorted_dyn_array_free(sorted_dyn_array* A);

typedef struct
{
    attr_id_t src;
    attr_id_t dest;
    attr_id_t eid;
} edge_t;

typedef struct
{
    dyn_array* neighbors;
    attr_id_t degree;
} vert_t;

/* The graph data structure*/
typedef struct
{ 
	/*** 
     The minimal graph repesentation consists of:
	 n        -- the number of vertices
	 m        -- the number of edges
	 endV     -- an array of size 2*m (m for directed graphs) that stores the 
                 destination ID of an edge <src->dest>.
	 numEdges -- an array of size n+1 that stores the degree 
                 (out-degree in case of directed graphs) and pointers to
                 the endV array. The degree of vertex i is given by 
                 numEdges[i+1]-numEdges[i], and the edges out of i are
                 stored in the contiguous block endV[numEdges[i] ..
                 numEdges[i+1]].
     Vertices are ordered from 0 in our internal representation
	 ***/
    long n;
    long m;
    attr_id_t* endV;
    attr_id_t* numEdges;
    
	int undirected;
    int zero_indexed;

	/***
     Other useful constructs that can be initialized when needed
     ***/
    edge_t* elist;
    edge_t* elist_aux;
    vert_t* vlist;
    vert_t* vlist_aux;

    /* For directed graphs, endBackV is used to store edges pointing into a
       vertex, and numBackEdges the in-degree */
    attr_id_t* endBackV;
    attr_id_t* numBackEdges;
    
    /* Vertex weights */
    int* int_weight_v;
    float* fl_weight_v;
    double* dbl_weight_v;
    long* l_weight_v;

    /* Edge weights */
    int weight_type;
    int* int_weight_e;
    float* fl_weight_e;
    double* dbl_weight_e;
    long* l_weight_e;

    double min_weight;
    double max_weight;

    /* Fine-grained locking support, currently using 
       OpenMP Mutex locks */
#ifdef _OPENMP
    omp_lock_t* vLock;
    omp_lock_t* eLock;
#endif

} graph_t;

typedef struct {
    attr_id_t* list;
    attr_id_t count;
    attr_id_t degree;
} plist;

#endif

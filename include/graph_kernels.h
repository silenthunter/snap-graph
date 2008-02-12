#ifndef _GRAPH_KERNELS_H
#define _GRAPH_KERNELS_H
#include "graph_defs.h"

long graph_diameter(graph_t* G);
void connected_components(graph_t* G, attr_id_t* component_num);
void strongly_connected_components(graph_t* G, attr_id_t* component_num);

void biconnected_components(graph_t* G, attr_id_t* component_num);
void find_articulation_points(graph_t* G, attr_id_t* component_num);
void biconnected_components_recursive(graph_t* G, attr_id_t v, attr_id_t p_v);
void art_points_recursive(graph_t* G, attr_id_t u);
void BiCC_stack_push(attr_id_t a, attr_id_t b);
void BiCC_stack_pop(attr_id_t *a, attr_id_t *b);

void BFS_parallel_frontier_expansion(graph_t* G, long src, long diameter);
void BFS_path_limited_search(graph_t* G, long src, long diameter);
void BFS_sequential(graph_t* G, long src, int* d);

#endif

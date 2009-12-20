#ifndef _GRAPH_PARTITIONING_H
#define _GRAPH_PARTITIONING_H
#include "graph_defs.h"
void modularity_spectral(graph_t *G, attr_id_t *membership, 
        attr_id_t *numCommunities, attr_id_t use_improvement);
void modularity_spectral_wo_klin(graph_t *G, attr_id_t *membership, 
        attr_id_t *numCommunities);
void computeModularityValue(graph_t *G, attr_id_t *membership, 
        attr_id_t numCommunities, double *modularity);
void modularity_betweenness(graph_t *g, attr_id_t *membership, 
        attr_id_t *num_communities, double *modularity, 
        double sampling_val);
void modularity_greedy_agglomerative(graph_t *g, char *alg_type, 
        attr_id_t *membership, attr_id_t *num_communities, double *modularity);
void single_seed_community_detection (graph_t *g, char *alg_type, attr_id_t seed, 
        attr_id_t *membership, attr_id_t *num_communities, double *modularity);
void seed_set_community_detection (graph_t *g, char *alg_type, attr_id_t* seeds, 
        int uniform, attr_id_t n_seeds, attr_id_t *membership, attr_id_t *num_communities,
		double *modularity);
long BFS_Community(graph_t *g, attr_id_t seed, long size, 
    attr_id_t* membership, attr_id_t *num_communities);
attr_id_t community_edges_external (graph_t *g, attr_id_t *membership, attr_id_t commId);
attr_id_t community_edges_internal (graph_t *g, attr_id_t *membership, attr_id_t commId);
#endif

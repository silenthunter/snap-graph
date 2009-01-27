#ifndef _GRAPH_PARTITIONING_H
#define _GRAPH_PARTITIONING_H

void modularity_spectral(graph_t *G, attr_id_t **membership, attr_id_t *numCommunities, attr_id_t use_improvement);
void computeModularityValue(graph_t *G, attr_id_t *membership, attr_id_t numCommunities, double *modularity);
void modularity_spectral_wo_klin(graph_t *G, attr_id_t **membership, attr_id_t *numCommunities);




typedef struct
{
        double *dist;   //stores the value
        int *vertex;    //map from the index of dist to vertex
        int *map;       //map from vertex to index of dist
        int *from;      // map from vertex to its parents in djkstra
        int size;
}heap;




#endif

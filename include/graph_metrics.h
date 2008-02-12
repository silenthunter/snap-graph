#ifndef _GRAPH_METRICS_H
#define _GRAPH_METRICS_H

void vertex_betweenness_centrality(graph_t* G, double* BC, long numSrcs);
void vertex_betweenness_centrality_simple(graph_t* G, double* BC, long numSrcs);
void vertex_betweenness_centrality_parBFS(graph_t* G, double* BC, long numSrcs);

void edge_betweenness_centrality(graph_t* G, double* BC, long numSrcs);
void edge_betweenness_centrality_simple(graph_t* G, double* BC, long numSrcs);
void edge_betweenness_centrality_parBFS(graph_t* G, double* BC, long numSrcs);

void vertex_closeness_centrality(graph_t* G, double* CC, long numSrcs);
void vertex_closeness_centrality_simple(graph_t* G, double* CC, long numSrcs);
void vertex_closeness_centrality_parBFS(graph_t* G, double* CC, long numSrcs);

void edge_closeness_centrality(graph_t* G, double* CC, long numSrcs);
void edge_closeness_centrality_simple(graph_t* G, double* CC, long numSrcs);
void edge_closeness_centrality_parBFS(graph_t* G, double* CC, long numSrcs);

#endif

#include "graph_partitioning.h"
#include "graph_kernels.h"

double comm_evaluate_modularity(graph_t* g, comm_list_t* comm_list, attr_id_t num_components) {

    attr_id_t u, v, j;
    attr_id_t n, m;
    attr_id_t cid;
    double m2_inv; 
    attr_id_t u_edge_start, u_edge_end, u_degree;
    double mod;

    n = g->n;
    m = g->m; 
    m2_inv = 1.0/(2 * g->m);
    mod = 0.0;
    
    for (j=0; j<num_components; j++) {
        comm_list[j].a = 0.0;
        comm_list[j].e = 0.0;
    }

    for(u=0; u<n; u++) {
        cid = g->cvl[u].comm_id;
        u_edge_start = g->cvl[u].num_edges;
        u_edge_end   = g->cvl[u+1].num_edges;
        u_degree     = u_edge_end - u_edge_start;

        for(j=u_edge_start; j<u_edge_end; j++) {
            v = g->cel[j].dest;
            if (g->cvl[v].comm_id == cid)
                comm_list[cid].e += 1.0;
        }
        comm_list[cid].a += (double) u_degree;
    }

    for (j=0; j<num_components; j++) {
        mod = comm_list[j].e - m2_inv * comm_list[j].a * comm_list[j].a; 
    }
    mod = mod * m2_inv;
    return mod;
}

void  remove_maxbc_edge(graph_t *g, comm_list_t *comm_list, 
        attr_id_t num_components, attr_id_t num_bc_runs, 
        edge_t* ebc_edge, attr_id_t* maxbc_component) {

    attr_id_t i, j;
    double mbc_val;
    attr_id_t mbc_eid;
    attr_id_t mbc_esrc, mbc_edest;
    attr_id_t mbc_component;
    attr_id_t i_start_edge, i_end_edge;
    attr_id_t n;

    n = g->n;

    /* find edge with max edge bc value */ 
    mbc_val = comm_list[0].mbc_val;
    for (i=1; i<num_components; i++) {
        if (comm_list[i].mbc_val > mbc_val) {     
            mbc_val  = comm_list[i].mbc_val;
            mbc_component = i;
            mbc_eid  = comm_list[i].mbc_eid;
            mbc_esrc = comm_list[i].mbc_esrc;     
        }
    }

    /* Mark the corresponding edge as deleted */
    i = mbc_esrc;
    i_start_edge = g->cvl[i].num_edges;
    i_end_edge = g->cvl[i+1].num_edges;
    for (j=i_start_edge; j<i_end_edge; j++) {
        if (g->cel[j].eid == mbc_eid) {
            mbc_edest = g->cel[j].dest; 
            g->cel[j].mask = num_bc_runs;
        } 
    }

    i = mbc_edest;
    i_start_edge = g->cvl[i].num_edges;
    i_end_edge = g->cvl[i+1].num_edges;
    for (j=i_start_edge; j<i_end_edge; j++) {
        if (g->cel[j].eid == mbc_eid) {
            g->cel[j].mask = num_bc_runs;
        } 
    }
    
    /* Update the maxbc edge in component */
    comm_list[mbc_component].mbc_val = -1;
    mbc_val = -1;
    for (i=0; i<n; i++) {
        if (g->cvl[i].comm_id != mbc_component)
            continue;
        i_start_edge = g->cvl[i].num_edges;
        i_end_edge   = g->cvl[i+1].num_edges;
        for (j=i_start_edge; j<i_end_edge; j++) {
            if (g->cel[j].mask == 0) {
                if (g->cel[j].cval > mbc_val) {    
                    mbc_val = g->cel[j].cval;
                    comm_list[mbc_component].mbc_val = mbc_val;
                    comm_list[mbc_component].mbc_eid = g->cel[j].eid;
                    comm_list[mbc_component].mbc_esrc = i;        
                }
            }
        }
    }

    *maxbc_component = mbc_component;

}
     
/* Community detection algorithm that optimizes modularity, and is based on 
   iterative removal of edges with high betweenness in the network */
void modularity_betweenness(graph_t *g, attr_id_t *membership, 
        attr_id_t *num_communities, double *modularity, 
        double sampling_val) {

    attr_id_t n, m;
    long i, j;    
    attr_id_t init_num_components, num_components, curr_component;
    double curr_modularity, prev_modularity;
    comm_list_t* comm_list;
    int new_comp;
    edge_t ebc_edge;
    int num_bc_runs;
    attr_id_t curr_component1, curr_component2, maxbc_component;
    int split;
    attr_id_t* ebc_eval_data1;
    double* ebc_eval_data2;

    n = g->n;
    m = g->m;
  
    /* Initialize the graph representation we will use in this algorithm */ 
    g->cvl = (c_vert_t *) calloc(n + 1, sizeof(c_vert_t));
    g->cel = (c_edge_t *) calloc(m, sizeof(c_edge_t));
    for (i=0; i<n; i++) {
        g->cvl[i].num_edges = g->numEdges[i];
        for (j=g->numEdges[i]; j<g->numEdges[i+1]; j++) {
            g->cel[j].dest = g->endV[j];
            g->cel[j].eid  = g->edge_id[j];
        }
    }

    g->cvl[n].num_edges = g->numEdges[n];
    /* Initially, all vertices belong to one giant community */
    fprintf(stderr, "Running connected components ...\n");

    /* Run connected components */
    num_components = aux_connected_components_init(g);
    fprintf(stderr, "The network has %d connected components.\n",
            num_components);

    init_num_components = num_components;
    
    /* We store the community splits to reconstruct the hierarchical 
       community dendrogram */
    /* Every community stores the ID of its parent, and also 
     * two variables to count the  */ 
    /* The max. number of communities is n, the number of vertices in the
     * network */
    comm_list = (comm_list_t *) calloc(n, sizeof(comm_list_t));
    curr_modularity = comm_evaluate_modularity(g, comm_list, num_components);
    prev_modularity = 0;
    num_bc_runs = 0;

    /* Initially run betweenness computation for all connected components */
    curr_component1 = curr_component2 = -1;
    split = 0;

    /* Pre-allocate memory for the edge BC computation in every iteration */

    /* start the iterative edge-betweenness based partitioning */
    while (1) {

        evaluate_edge_centrality_bcpart(g, ebc_eval_data1, ebc_eval_data2, comm_list, num_components, curr_component1, curr_component2);
        remove_maxbc_edge(g, comm_list, num_components, num_bc_runs, &ebc_edge, 
                &maxbc_component);
        num_bc_runs++;

        split = aux_connected_components_update(g, num_components,
                maxbc_component);
        
        if (split) {
            curr_component1 = maxbc_component;
            curr_component2 = num_components;
            comm_list[num_components].p = maxbc_component;
            num_components++;
            prev_modularity = curr_modularity;
            curr_modularity = comm_evaluate_modularity(g, comm_list, num_components);
            if (curr_modularity < prev_modularity)
                break;
        }

    } 

}

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

/**
 * \function igraph_community_edge_betweenness
 * \brief Community findinf based on edge betweenness
 * 
 * Community structure detection based on the betweenness of the edges
 * in the network. The algorithm was invented by M. Girvan and
 * M. Newman, see: M. Girvan and M. E. J. Newman: Community structure in
 * social and biological networks, Proc. Nat. Acad. Sci. USA 99, 7821-7826
 * (2002).
 * 
 * </para><para>
 * The idea is that the betweenness of the edges connecting two
 * communities is typically high, as many of the shortest paths
 * between nodes in separate communities go through them. So we
 * gradually remove the edge with highest betweenness from the
 * network, and recalculate edge betweenness after every removal. 
 * This way sooner or later the network falls off to two components,
 * then after a while one of these components falls off to two smaller 
 * components, etc. until all edges are removed. This is a divisive
 * hieararchical approach, the result is a dendrogram.
 * \param graph The input graph.
 * \param result Pointer to an initialized vector, the result will be
 *     stored here, the ids of the removed edges in the order of their 
 *     removal. It will be resized as needed.
 * \param edge_betweenness Pointer to an initialized vector or
 *     NULL. In the former case the edge betweenness of the removed
 *     edge is stored here. The vector will be resized as needed.
 * \param merges Pointer to an initialized matrix or NULL. If not NULL
 *     then merges performed by the algorithm are stored here. Even if
 *     this is a divisive algorithm, we can replay it backwards and
 *     note which two clusters were merged. Clusters are numbered from
 *     zero, see the \p merges argument of \ref
 *     igraph_community_walktrap() for details. The matrix will be
 *     resized as needed.
 * \param bridges Pointer to an initialized vector of NULL. If not
 *     NULL then all edge removals which separated the network into
 *     more components are marked here.
 * \param directed Logical constant, whether to calculate directed
 *    betweenness (ie. directed paths) for directed graphs. It is
 *    ignored for undirected graphs.
 * \return Error code.
 * 
 * \sa \ref igraph_community_eb_get_merges(), \ref
 * igraph_community_spinglass(), \ref igraph_community_walktrap().
 * 
 * Time complexity: O(|V|^3), as the betweenness calculation requires
 * O(|V|^2) and we do it |V|-1 times.
 */
#if 0  
int igraph_community_edge_betweenness(const igraph_t *graph, 
				      igraph_vector_t *result,
				      igraph_vector_t *edge_betweenness,
				      igraph_matrix_t *merges,
				      igraph_vector_t *bridges,
				      igraph_bool_t directed) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance, *nrgeo;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source, i, e;
  
  igraph_adjedgelist_t elist_out, elist_in;
  igraph_adjedgelist_t *elist_out_p, *elist_in_p;
  igraph_vector_t *neip;
  long int neino;
  igraph_integer_t modein, modeout;
  igraph_vector_t eb;
  long int maxedge, pos;
  igraph_integer_t from, to;

  char *passive;

  directed=directed && igraph_is_directed(graph);
  if (directed) {
    modeout=IGRAPH_OUT;
    modeout=IGRAPH_IN;
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_out);
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_in);
    elist_out_p=&elist_out;
    elist_in_p=&elist_in;
  } else {
    modeout=modein=IGRAPH_ALL;
    IGRAPH_CHECK(igraph_adjedgelist_init(graph, &elist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjedgelist_destroy, &elist_out);
    elist_out_p=elist_in_p=&elist_out;
  }
  
  distance=igraph_Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, distance);
  nrgeo=igraph_Calloc(no_of_nodes, long int);
  if (nrgeo==0) {
    IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, nrgeo);
  tmpscore=igraph_Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, tmpscore);

  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_CHECK(igraph_stack_init(&stack, no_of_nodes));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
  
  IGRAPH_CHECK(igraph_vector_resize(result, no_of_edges));
  if (edge_betweenness) {
    IGRAPH_CHECK(igraph_vector_resize(edge_betweenness, no_of_edges));
    VECTOR(*edge_betweenness)[no_of_edges-1]=0;
  }

  IGRAPH_VECTOR_INIT_FINALLY(&eb, no_of_edges);
  
  passive=igraph_Calloc(no_of_edges, char);
  if (!passive) {
    IGRAPH_ERROR("edge betweenness community structure failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, passive);

  for (e=0; e<no_of_edges; e++) {
    
    igraph_vector_null(&eb);

    for (source=0; source<no_of_nodes; source++) {

      /* This will contain the edge betweenness in the current step */
      IGRAPH_ALLOW_INTERRUPTION();

      memset(distance, 0, no_of_nodes*sizeof(long int));
      memset(nrgeo, 0, no_of_nodes*sizeof(long int));
      memset(tmpscore, 0, no_of_nodes*sizeof(double));
      igraph_stack_clear(&stack); /* it should be empty anyway... */
      
      IGRAPH_CHECK(igraph_dqueue_push(&q, source));
      
      nrgeo[source]=1;
      distance[source]=0;
      
      while (!igraph_dqueue_empty(&q)) {
	long int actnode=igraph_dqueue_pop(&q);
	
	neip=igraph_adjedgelist_get(elist_out_p, actnode);
	neino=igraph_vector_size(neip);
	for (i=0; i<neino; i++) {
	  igraph_integer_t edge=VECTOR(*neip)[i], from, to;
	  long int neighbor;
	  igraph_edge(graph, edge, &from, &to);
	  neighbor = actnode!=from ? from : to;
	  if (nrgeo[neighbor] != 0) {
	    /* we've already seen this node, another shortest path? */
	    if (distance[neighbor]==distance[actnode]+1) {
	      nrgeo[neighbor]+=nrgeo[actnode];
	    }
	  } else {
	    /* we haven't seen this node yet */
	    nrgeo[neighbor]+=nrgeo[actnode];
	    distance[neighbor]=distance[actnode]+1;
	    IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	    IGRAPH_CHECK(igraph_stack_push(&stack, neighbor));
	  }
	}
      } /* while !igraph_dqueue_empty */
      
      /* Ok, we've the distance of each node and also the number of
	 shortest paths to them. Now we do an inverse search, starting
	 with the farthest nodes. */
      while (!igraph_stack_empty(&stack)) {
	long int actnode=igraph_stack_pop(&stack);
	if (distance[actnode]<1) { continue; } /* skip source node */
	
	/* set the temporary score of the friends */
	neip=igraph_adjedgelist_get(elist_in_p, actnode);
	neino=igraph_vector_size(neip);
	for (i=0; i<neino; i++) {
	  igraph_integer_t from, to;
	  long int neighbor;
	  long int edgeno=VECTOR(*neip)[i];
	  igraph_edge(graph, edgeno, &from, &to);
	  neighbor= actnode != from ? from : to;
	  if (distance[neighbor]==distance[actnode]-1 &&
	      nrgeo[neighbor] != 0) {
	    tmpscore[neighbor] +=
	      (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	    VECTOR(eb)[edgeno] +=
	      (tmpscore[actnode]+1)*nrgeo[neighbor]/nrgeo[actnode];
	  }
	}
      }
      /* Ok, we've the scores for this source */
    } /* for source <= no_of_nodes */
    
    /* Now look for the smallest edge betweenness */
    /* and eliminate that edge from the network */
    maxedge=igraph_i_vector_which_max_not_null(&eb, passive);
    VECTOR(*result)[e]=maxedge;
    if (edge_betweenness) {
      VECTOR(*edge_betweenness)[e]=VECTOR(eb)[maxedge];
      if (!directed) { 
	VECTOR(*edge_betweenness)[e] /= 2.0;
      }
    }
    passive[maxedge]=1;
    igraph_edge(graph, maxedge, &from, &to);

    neip=igraph_adjedgelist_get(elist_in_p, to);
    neino=igraph_vector_size(neip);
    igraph_vector_search(neip, 0, maxedge, &pos);
    VECTOR(*neip)[pos]=VECTOR(*neip)[neino-1];
    igraph_vector_pop_back(neip);
    
    neip=igraph_adjedgelist_get(elist_out_p, from);
    neino=igraph_vector_size(neip);
    igraph_vector_search(neip, 0, maxedge, &pos);
    VECTOR(*neip)[pos]=VECTOR(*neip)[neino-1];
    igraph_vector_pop_back(neip);
  }

  igraph_free(passive);
  igraph_vector_destroy(&eb);
  igraph_stack_destroy(&stack);
  igraph_dqueue_destroy(&q);
  igraph_free(tmpscore);
  igraph_free(nrgeo);
  igraph_free(distance);
  IGRAPH_FINALLY_CLEAN(7);
  
  if (directed) {
    igraph_adjedgelist_destroy(&elist_out);
    igraph_adjedgelist_destroy(&elist_in);
    IGRAPH_FINALLY_CLEAN(2);
  } else {
    igraph_adjedgelist_destroy(&elist_out);
    IGRAPH_FINALLY_CLEAN(1);
  }

  if (merges || bridges) {
    IGRAPH_CHECK(igraph_community_eb_get_merges(graph, result, merges, bridges));
  }
  
  return 0;
}


/**
 * \function igraph_community_to_membership 
 * \brief Create membership vector from community structure dendrogram
 * 
 * This function creates a membership vector from a community
 * structure dendrogram. A membership vector contains for each vertex
 *
 * 
 */
#endif

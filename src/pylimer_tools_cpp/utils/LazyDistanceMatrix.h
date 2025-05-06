#pragma once

extern "C"
{
#include <igraph/igraph.h>
}
#include <cassert>

typedef struct igraph_lazy_distance_matrix_state_t
{
  /* Number of vertices in the graph */
  igraph_integer_t N;

  /* The adjacency list */
  igraph_lazy_adjlist_t adjlist;

  /* Whether the graph is directed */
  igraph_bool_t directed;

  /* The path distance cut-off */
  igraph_integer_t max_path_length;

  /* The resulting path lengths */
  igraph_vector_int_list_t path_lengths;

  /* A vector containing the vertices that have been visited in the current
   * search iteration */
  igraph_vector_int_t visiting_hosts_in_current_iteration;

  /* A vector indicating for each vertex whether it has been visited in the
   * current search iteration */
  igraph_vector_bool_t is_visiting_in_current_interation;

  /*  */
  igraph_dqueue_int_t to_visit_queue;
  igraph_dqueue_int_t depth_queue;
} igraph_lazy_distance_matrix_state_t;

static void
igraph_lazy_distance_matrix_forget_source(
  const igraph_lazy_distance_matrix_state_t* state,
  const igraph_integer_t source)
{
  igraph_vector_int_t* path_lengths =
    igraph_vector_int_list_get_ptr(&state->path_lengths, source);
  igraph_vector_int_resize(path_lengths, 0);
  igraph_vector_int_resize_min(path_lengths);
}

/**
 * @brief Retrieves the path length between two vertices in a graph
 *
 * This function calculates the shortest path length between source and target
 * vertices using breadth-first search. It can optionally store all computed
 * path lengths.
 *
 * @param state The lazy distance matrix state containing graph information
 * @param source The source vertex ID
 * @param target The target vertex ID
 * @param store_path Whether to store all path lengths (default: false)
 * @return The length of the shortest path between source and target,
 * or IGRAPH_INFINITY if no path exists
 */
static igraph_integer_t
igraph_lazy_distance_matrix_path_length(
  igraph_lazy_distance_matrix_state_t* state,
  const igraph_integer_t source,
  const igraph_integer_t target,
  const igraph_bool_t store_path = false)
{
  // we only need one direction in the case of undirected graphs
  // make sure we query always the same matrix element
  if (!state->directed && source > target) {
    return igraph_lazy_distance_matrix_path_length(state, target, source);
  }
  // check whether we have already computed the path length
  igraph_vector_int_t* path_lengths =
    igraph_vector_int_list_get_ptr(&state->path_lengths, source);
  if (store_path) {
    if (igraph_vector_int_size(path_lengths) != state->N) {
      igraph_vector_int_resize(path_lengths, state->N);
      igraph_vector_int_null(path_lengths);
    }
    if (igraph_vector_int_get(path_lengths, target) > 0) {
      return igraph_vector_int_get(path_lengths, target);
    }
  }
  // otherwise, compute it by bfs
  igraph_dqueue_int_push(&state->to_visit_queue, source);
  igraph_dqueue_int_push(&state->depth_queue, 0);
  igraph_integer_t result = IGRAPH_INFINITY;
  // remember having visited the start
  igraph_vector_int_push_back(&state->visiting_hosts_in_current_iteration,
                              source);
  VECTOR(state->is_visiting_in_current_interation)
  [source] = true;
  while (igraph_dqueue_int_size(&state->to_visit_queue) > 0) {
    // fetch the next vertex and its depth
    const igraph_integer_t current_vertex =
      igraph_dqueue_int_pop(&state->to_visit_queue);
    const igraph_integer_t current_depth =
      igraph_dqueue_int_pop(&state->depth_queue);

    // iterate the neighbors to either walk deeper, or find target vertex
    const igraph_vector_int_t* neighbors =
      igraph_lazy_adjlist_get(&state->adjlist, current_vertex);
    const igraph_integer_t next_depth = current_depth + 1;
    for (igraph_integer_t neighbour_idx = 0;
         neighbour_idx < igraph_vector_int_size(neighbors);
         ++neighbour_idx) {
      const igraph_integer_t neighbour_vertex =
        igraph_vector_int_get(neighbors, neighbour_idx);
      // store the path length in any case, since we have already computed it
      if (store_path) {
        const igraph_integer_t dist =
          igraph_vector_int_get(path_lengths, neighbour_vertex);
        assert(dist == 0 || dist == next_depth);
        igraph_vector_int_set(path_lengths, neighbour_vertex, next_depth);
      }
      if (neighbour_vertex == target) {
        result = next_depth;
        break;
      }
      if (!VECTOR(state->is_visiting_in_current_interation)[neighbour_vertex] &&
          (next_depth) < state->max_path_length) {
        igraph_dqueue_int_push(&state->to_visit_queue, neighbour_vertex);
        igraph_dqueue_int_push(&state->depth_queue, next_depth);
        // remember having visited this vertex, in order not to visit it again
        igraph_vector_int_push_back(&state->visiting_hosts_in_current_iteration,
                                    current_vertex);
        VECTOR(state->is_visiting_in_current_interation)
        [current_vertex] = true;
      }
    }
    assert(igraph_dqueue_int_size(&state->to_visit_queue) ==
           igraph_dqueue_int_size(&state->depth_queue));
    // C does not have a `break 2`, nor should one use a `goto` statement,
    // so we break again here if we found the target vertex
    if (result != IGRAPH_INFINITY) {
      break;
    }
  }
  // reset visited flags for next iteration
  for (igraph_integer_t visited_idx = 0;
       visited_idx <
       igraph_vector_int_size(&state->visiting_hosts_in_current_iteration);
       ++visited_idx) {
    VECTOR(state->is_visiting_in_current_interation)
    [VECTOR(state->visiting_hosts_in_current_iteration)[visited_idx]] = false;
  }
  igraph_vector_int_clear(&state->visiting_hosts_in_current_iteration);
  igraph_dqueue_int_clear(&state->to_visit_queue);
  igraph_dqueue_int_clear(&state->depth_queue);

  // finally, return results
  return result;
}

static igraph_error_t
igraph_lazy_distance_matrix_state_init(
  igraph_lazy_distance_matrix_state_t* state,
  const igraph_t* graph,
  const igraph_neimode_t mode,
  const igraph_integer_t max_path_length)
{
  state->directed = igraph_is_directed(graph) && mode != IGRAPH_ALL;
  state->max_path_length = max_path_length;
  state->N = igraph_vcount(graph);
  IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&state->path_lengths, state->N);
  IGRAPH_VECTOR_INT_INIT_FINALLY(&state->visiting_hosts_in_current_iteration,
                                 state->N);
  igraph_vector_int_clear(&state->visiting_hosts_in_current_iteration);
  IGRAPH_VECTOR_BOOL_INIT_FINALLY(&state->is_visiting_in_current_interation,
                                  state->N);
  IGRAPH_DQUEUE_INT_INIT_FINALLY(&state->to_visit_queue, state->N);
  IGRAPH_DQUEUE_INT_INIT_FINALLY(&state->depth_queue, state->N);
  IGRAPH_CHECK(igraph_lazy_adjlist_init(
    graph, &state->adjlist, mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));

  IGRAPH_FINALLY_CLEAN(5);
  return IGRAPH_SUCCESS;
}

static void
igraph_lazy_distance_matrix_state_destroy(
  igraph_lazy_distance_matrix_state_t* state)
{
  igraph_lazy_adjlist_destroy(&state->adjlist);
  igraph_vector_int_list_destroy(&state->path_lengths);
  igraph_vector_int_destroy(&state->visiting_hosts_in_current_iteration);
  igraph_vector_bool_destroy(&state->is_visiting_in_current_interation);
  igraph_dqueue_int_destroy(&state->to_visit_queue);
  igraph_dqueue_int_destroy(&state->depth_queue);
}

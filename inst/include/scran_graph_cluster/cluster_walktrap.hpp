#ifndef SCRAN_CLUSTER_WALKTRAP_HPP
#define SCRAN_CLUSTER_WALKTRAP_HPP

#include <vector>
#include <algorithm>

#include "raiigraph/raiigraph.hpp"
#include "igraph.h"

/**
 * @file cluster_walktrap.hpp
 * @brief Wrapper around **igraph**'s Walktrap community detection algorithm.
 */

namespace scran_graph_cluster {

/**
 * @brief Options for `cluster_walktrap()`.
 */
struct ClusterWalktrapOptions {
    /**
     * Number of steps of the random walk.
     * The default is based on the example in the **igraph** documentation.
     */
    int steps = 4;

    /**
     * Whether to report the merge steps in `Results::merges`.
     */
    bool report_merges = true;

    /**
     * Whether to report the modularity after each merge step in `Results::modularity`.
     */
    bool report_modularity = true;
};

/**
 * @brief Result of `cluster_walktrap()`.
 */
struct ClusterWalktrapResults {
    /** 
     * Output status.
     * A value of zero indicates that the algorithm completed successfully.
     */
    int status = 0;
    
    /**
     * Vector of length equal to the number of cells, containing 0-indexed cluster identities.
     */
    raiigraph::IntegerVector membership;

    /**
     * Matrix of merge steps.
     * Each row corresponds to a successive merge step, while the two columns contain the identities of the two clusters being merged.
     * Note that cluster IDs here are not the same as those in `membership` - 
     * see [the documentation](https://igraph.org/c/doc/igraph-Community.html#igraph_community_walktrap) for more details.
     * This should only be used if `ClusterWalktrapOptions::report_merges = true`.
     */
    raiigraph::IntegerMatrix merges;

    /**
     * Vector of length equal to the number of rows in `merges` plus 1, containing the modularity score before and after each merge step.
     * The maximum value is the modularity corresponding to the clustering in `membership`.
     * This should only be used if `ClusterWalktrapOptions::report_modularity = true`.
     */
    raiigraph::RealVector modularity;
};

/**
 * Run the Walktrap community detection algorithm on a pre-constructed graph to obtain communities of highly inter-connected nodes.
 * See [here](https://igraph.org/c/doc/igraph-Community.html#igraph_community_walktrap) for more details on the Walktrap algorithm. 
 * 
 * @param graph An existing graph.
 * @param weights Pointer to an array of weights of length equal to the number of edges in `graph`. 
 * This should be in the same order as the edge list in `graph`.
 * Alternatively `NULL`, if the graph is unweighted.
 * @param options Further options.
 * @param[out] output On output, this is filtered with the clustering results.
 * The input value is ignored, so this object can be re-used across multiple calls to `cluster_walktrap()`.
 */
inline void cluster_walktrap(const igraph_t* graph, const igraph_vector_t* weights, const ClusterWalktrapOptions& options, ClusterWalktrapResults& output) {
    auto membership = output.membership.get();
    auto modularity = (options.report_modularity ? output.modularity.get() : NULL);
    auto merges = (options.report_merges ? output.merges.get() : NULL);
    output.status = igraph_community_walktrap(graph, weights, options.steps, merges, modularity, membership);
}

/**
 * Overload of `cluster_walktrap()` that accepts C++ containers instead of the raw **igraph** pointers.
 *
 * @param graph An existing graph.
 * @param weights Vector of weights of length equal to the number of edges in `graph`. 
 * This should be in the same order as the edge list in `graph`.
 * @param options Further options.
 *
 * @return Clustering results for the nodes of the graph.
 */
inline ClusterWalktrapResults cluster_walktrap(const raiigraph::Graph& graph, const std::vector<igraph_real_t>& weights, const ClusterWalktrapOptions& options) {
    // No need to free this, as it's just a view.
    igraph_vector_t weight_view;
    igraph_vector_view(&weight_view, weights.data(), weights.size());

    ClusterWalktrapResults output;
    cluster_walktrap(graph.get(), &weight_view, options, output);
    return output;
}

}

#endif

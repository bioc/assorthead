#ifndef SCRAN_GRAPH_CLUSTER_CLUSTER_LEIDEN_HPP
#define SCRAN_GRAPH_CLUSTER_CLUSTER_LEIDEN_HPP

#include <vector>
#include <algorithm>

#include "raiigraph/raiigraph.hpp"
#include "igraph.h"

/**
 * @file cluster_leiden.hpp
 * @brief Wrapper around **igraph**'s Leiden community detection algorithm.
 */

namespace scran_graph_cluster {

/**
 * @brief Options for `cluster_leiden()`.
 */
struct ClusterLeidenOptions {
    /**
     * Resolution of the clustering.
     * Larger values result in more fine-grained communities.
     * The default is based on `?cluster_leiden` in the **igraph** R package.
     */
    double resolution = 1;

    /**
     * Level of randomness used during refinement.
     * The default is based on `?cluster_leiden` in the **igraph** R package.
     */
    double beta = 0.01;

    /**
     * Number of iterations of the Leiden algorithm.
     * More iterations can improve separation at the cost of computational time.
     * The default is based on `?cluster_leiden` in the **igraph** R package.
     */
    int iterations = 2;

    /**
     * Whether to optimize the modularity instead of the Constant Potts Model.
     * The two are closely related but the magnitude of the resolution is different.
     * The default is based on `?cluster_leiden` in the **igraph** R package.
     */
    bool modularity = false;

    /**
     * Whether to report the quality of the clustering in `Results::quality`.
     */
    bool report_quality = true;

    /**
     * Seed for the **igraph** random number generator.
     */
    int seed = 42;
};

/**
 * @brief Result of `cluster_leiden()`.
 */
struct ClusterLeidenResults {
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
     * Quality of the clustering, closely related to the modularity.
     * This is only defined if `ClusterLeidenOptions::report_quality = true`.
     */
    igraph_real_t quality = 0;
};

/**
 * Run the Leiden community detection algorithm on a pre-constructed graph to obtain communities of highly inter-connected nodes.
 * See [here](https://igraph.org/c/doc/igraph-Community.html#igraph_community_leiden) for more details on the Leiden algorithm. 
 *
 * @param graph An existing graph.
 * @param weights Pointer to an array of weights of length equal to the number of edges in `graph`. 
 * This should be in the same order as the edge list in `graph`.
 * Alternatively `NULL`, if the graph is unweighted.
 * @param options Further options.
 * @param[out] output On output, this is filtered with the clustering results.
 * The input value is ignored, so this object can be re-used across multiple calls to `cluster_leiden()`.
 */
inline void cluster_leiden(const igraph_t* graph, const igraph_vector_t* weights, const ClusterLeidenOptions& options, ClusterLeidenResults& output) {
    auto membership = output.membership.get();
    auto quality = (options.report_quality ? &(output.quality) : NULL);

    raiigraph::RNGScope rngs(options.seed);

    if (!options.modularity) {
        output.status = igraph_community_leiden(
            graph, 
            weights,
            NULL,
            options.resolution, 
            options.beta,
            false, 
            options.iterations, 
            membership, 
            NULL,
            quality
        );

    } else {
        // More-or-less translated from igraph::cluster_leiden in the R package.
        raiigraph::RealVector strengths(igraph_vcount(graph));
        igraph_strength(graph, strengths, igraph_vss_all(), IGRAPH_ALL, 1, weights);

        double total_weights = igraph_vector_sum(strengths);
        double mod_resolution = options.resolution / total_weights;

        output.status = igraph_community_leiden(
            graph, 
            weights, 
            strengths, 
            mod_resolution, 
            options.beta, 
            false, 
            options.iterations, 
            membership, 
            NULL, 
            quality
        );
    }
}

/**
 * Overload of `cluster_leiden()` that accepts C++ containers instead of the raw **igraph** pointers.
 *
 * @param graph An existing graph.
 * @param weights Vector of weights of length equal to the number of edges in `graph`. 
 * This should be in the same order as the edge list in `graph`.
 * @param options Further options.
 *
 * @return Clustering results for the nodes of the graph.
 */
inline ClusterLeidenResults cluster_leiden(const raiigraph::Graph& graph, const std::vector<igraph_real_t>& weights, const ClusterLeidenOptions& options) {
    // No need to free this, as it's just a view.
    igraph_vector_t weight_view;
    igraph_vector_view(&weight_view, weights.data(), weights.size());

    ClusterLeidenResults output;
    cluster_leiden(graph.get(), &weight_view, options, output);
    return output;
}

}

#endif

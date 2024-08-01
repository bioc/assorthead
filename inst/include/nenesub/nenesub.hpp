#ifndef NENESUB_HPP
#define NENESUB_HPP

#include <vector>
#include <queue>
#include <cstdint>
#include <algorithm>
#include "knncolle/knncolle.hpp"

/**
 * @file nenesub.hpp
 * @brief Nearest-neighbors subsampling.
 */

/**
 * @namespace nenesub
 * @brief Nearest-neighbors subsampling.
 */
namespace nenesub {

/**
 * @brief Options for `compute()`.
 */
struct Options {
    /**
     * The number of nearest neighbors to use, i.e., \f$k\f$. 
     * Larger values decrease the subsampling rate, i.e., fewer observations are selected.
     * Only relevant for the `compute()` overloads without pre-computed neighbors.
     */
    int num_neighbors = 20;

    /**
     * The minimum number of remaining neighbors that an observation must have in order to be selected.
     * Larger values decrease the subsampling rate, i.e., fewer observations are selected.
     * This should be less than or equal to `Options::num_neighbors`.
     */
    int min_remaining = 10;

    /**
     * The number of threads to use.
     * Only relevant for the `compute()` overloads without pre-computed neighbors.
     */
    int num_threads = 10;
};

/**
 * This function generates a deterministic subsampling of a dataset based on nearest neighbors.
 * We first identify the \f$k\f$-nearest neighbors of each observation and use that to define its local neighborhood.
 * We select an observation for subsampling if it:
 *
 * - Does not belong in the local neighborhood of any previously selected observation.
 * - Has the most neighbors that are not selected or in the local neighborhoods of previously selected observations.
 *   Ties are broken using the smallest distance to each observation's \f$k\f$-th neighbor (i.e., the densest region of space).
 * - Has at least `Options::min_remaining` neighbors that are not selected or in the local neighborhoods of any other selected observation.
 *
 * We repeat this process until there are no more observations that satisfy these requirements. 
 * Each selected observation serves as a representative for up to \f$k\f$ of its nearest neighbors.
 * As such, the rate of subsampling is roughly proportional to the chocie of \f$k\f$, e.g., \f$k = 20\f$ suggests that every 20th observation will be selected on average.
 *
 * The **nenesub** approach ensures that the subsampled points are well-distributed across the dataset.
 * Low-frequency subpopulations will always have at least a few representatives if they are sufficiently distant from other subpopulations.
 * In contrast, random sampling does not provide strong guarantees for capture of a rare subpopulation.
 * We also preserve the relative density across the dataset as more representatives will be generated from high-density regions. 
 * This simplifies the interpretation of analysis results generated from the subsetted dataset.
 * 
 * @tparam Index_ Integer type for the observation indices.
 * @tparam GetNeighbors_ Function that accepts an `Index_` index and returns a (const reference to a) container-like object.
 * The container should be support the `[]` operator and have a `size()` method.
 * @tparam GetIndex_ Function that accepts a (const reference to a) container of the type returned by `GetNeighbors_` and an `Index_` into that container, and returns `Index_`.
 * @tparam GetNeighbors_ Function that accepts an `Index_` index and returns a distance value, typically floating-point.
 *
 * @param num_obs Number of observations in the dataset.
 * @param get_neighbors Function that accepts an integer observation index in `[0, num_obs)` and returns a container of that observation's neighbors.
 * Each element of the container specifies the index of a neighboring observation.
 * It is generally expected that the returned containers have the same size for all indices.
 * @param get_index Function to return the index of each neighbor, given the container returned by `get_neighbors` and an index into that container.
 * @param get_max_distance Function that accepts an integer observation index in `[0, num_obs)` and returns the distance from that observation to its furthest neighbor.
 * @param options Further options. 
 * Note that `Options::num_neighbors` and `Options::num_threads` are ignored here.
 * @param[out] selected On output, the indices of the observations that were subsampled.
 * These are sorted in ascending order.
 */
template<typename Index_, class GetNeighbors_, class GetIndex_, class GetMaxDistance_>
void compute(Index_ num_obs, GetNeighbors_ get_neighbors, GetIndex_ get_index, GetMaxDistance_ get_max_distance, const Options& options, std::vector<Index_>& selected) {
    typedef decltype(get_max_distance(0)) Distance_;
    struct Payload {
        Payload(Index_ identity, Index_ remaining, Distance_ max_distance) : remaining(remaining), identity(identity), max_distance(max_distance) {}
        Index_ remaining;
        Index_ identity;
        Distance_ max_distance;
    };

    auto cmp = [](const Payload& left, const Payload& right) -> bool {
        if (left.remaining == right.remaining) {
            if (left.max_distance == right.max_distance) {
                return left.identity > right.identity; // smallest identities show up first.
            }
            return left.max_distance > right.max_distance; // smallest distances show up first.
        }
        return left.remaining < right.remaining; // largest remaining show up first.
    };
    std::priority_queue<Payload, std::vector<Payload>, decltype(cmp)> store(
        cmp,
        [&]{
            std::vector<Payload> container;
            container.reserve(num_obs);
            return container;
        }()
    );

    std::vector<std::vector<Index_> > reverse_map(num_obs);
    std::vector<Index_> remaining(num_obs);
    for (Index_ c = 0; c < num_obs; ++c) {
        const auto& neighbors = get_neighbors(c);
        Index_ nneighbors = neighbors.size();

        if (nneighbors) { // protect get_max_distance just in case there are no neighbors.
            store.emplace(c, nneighbors, get_max_distance(c));
            for (Index_ n = 0; n < nneighbors; ++n) {
                reverse_map[get_index(neighbors, n)].push_back(c);
            }
            remaining[c] = nneighbors;
        }
    }

    selected.clear();
    std::vector<uint8_t> tainted(num_obs);
    while (!store.empty()) {
        auto payload = store.top();
        store.pop();
        if (tainted[payload.identity]) {
            continue;
        }

        const auto& neighbors = get_neighbors(payload.identity);
        Index_ new_remaining = remaining[payload.identity];

        if (new_remaining >= options.min_remaining) {
            payload.remaining = new_remaining;
            if (!store.empty() && cmp(payload, store.top())) {
                store.push(payload);
            } else {
                selected.push_back(payload.identity);
                tainted[payload.identity] = 1;
                for (auto x : reverse_map[payload.identity]) {
                    --remaining[x];
                }

                Index_ nneighbors = neighbors.size();
                for (Index_ n = 0; n < nneighbors; ++n) {
                    auto current = get_index(neighbors, n);
                    tainted[current] = 1;
                    for (auto x : reverse_map[current]) {
                        --remaining[x];
                    }
                }
            }
        }
    }

    std::sort(selected.begin(), selected.end());
}

/**
 * Overload to enable convenient usage with pre-computed neighbors from **knncolle**.
 *
 * @tparam Index_ Integer type for the neighbor indices.
 * @tparam Distance_ Floating-point type for the distances.
 *
 * @param neighbors Vector of nearest-neighbor search results for each observation.
 * Each entry is a pair containing a vector of neighbor indices and a vector of distances to those neighbors.
 * Neighbors should be sorted by increasing distance.
 * The same number of neighbors should be present for each observation.
 * @param options Further options. 
 * Note that `Options::num_neighbors` and `Options::num_threads` are ignored here.
 *
 * @return A sorted vector of the indices of the subsampled observations.
 */
template<typename Index_, typename Distance_>
std::vector<Index_> compute(const knncolle::NeighborList<Index_, Distance_>& neighbors, const Options& options) {
    std::vector<Index_> output;
    compute(
        static_cast<Index_>(neighbors.size()),
        [&](size_t i) -> const std::vector<Index_>& { return neighbors[i].first; }, 
        [](const std::vector<Index_>& x, Index_ n) -> Index_ { return x[n]; }, 
        [&](size_t i) -> Distance_ { return neighbors[i].second.back(); }, 
        options,
        output
    );
    return output;
}

/**
 * Overload to enable convenient usage with a prebuilt nearest-neighbor search index from **knncolle**.
 *
 * @tparam Dim_ Integer type for the dimension index.
 * @tparam Index_ Integer type for the observation index.
 * @tparam Float_ Floating-point type for the distances.
 *
 * @param[in] prebuilt A prebuilt nearest-neighbor search index on the observations of interest.
 * @param options Further options.
 *
 * @return A sorted vector of the indices of the subsampled observations.
 */
template<typename Dim_, typename Index_, typename Float_>
std::vector<Index_> compute(const knncolle::Prebuilt<Dim_, Index_, Float_>& prebuilt, const Options& options) {
    Index_ nobs = prebuilt.num_observations();
    std::vector<std::vector<Index_> > nn_indices(nobs);
    std::vector<Float_> max_distance(nobs);
    int k = options.num_neighbors;

#ifndef KNNCOLLE_CUSTOM_PARALLEL
#ifdef _OPENMP
    #pragma omp parallel num_threads(options.num_threads)
    {
    auto sptr = prebuilt.initialize();
    std::vector<Float_> nn_distances;
    #pragma omp for
    for (Index_ i = 0; i < nobs; ++i) {
#else
    auto sptr = prebuilt.initialize();
    std::vector<Float_> nn_distances;
    for (Index_ i = 0; i < nobs; ++i) {
#endif
#else
    KNNCOLLE_CUSTOM_PARALLEL(nobs, options.num_threads, [&](Index_ start, Index_ length) -> void {
    auto sptr = prebuilt.initialize();
    std::vector<Float_> nn_distances;
    for (Index_ i = start, end = start + length; i < end; ++i) {
#endif        

        sptr->search(i, k, &(nn_indices[i]), &nn_distances);
        max_distance[i] = (k ? 0 : nn_distances.back());

#ifndef KNNCOLLE_CUSTOM_PARALLEL    
#ifdef _OPENMP
    }
    }
#else
    }
#endif
#else
    }
    });
#endif

    std::vector<Index_> output;
    compute(
        nobs,
        [&](size_t i) -> const std::vector<Index_>& { return nn_indices[i]; }, 
        [](const std::vector<Index_>& x, Index_ n) -> Index_ { return x[n]; }, 
        [&](size_t i) -> Float_ { return max_distance[i]; },
        options,
        output
    );
    return output;
}

/**
 * Overload to enable convenient usage with a column-major array of coordinates for each observation.
 *
 * @tparam Dim_ Integer type for the dimension index.
 * @tparam Index_ Integer type for the observation index.
 * @tparam Value_ Numeric type for the input data.
 * @tparam Float_ Floating-point type for the distances.
 *
 * @param num_dims Number of dimensions for the observation coordinates.
 * @param num_obs Number of observations in the dataset.
 * @param[in] data Pointer to a `num_dims`-by-`num_observations` column-major array of observation coordinates where rows are dimensions and columns are observations.
 * @param knn_method Specification of the nearest-neighbor search algorithm, e.g., `knncolle::VptreeBuilder`, `knncolle::KmknnBuilder`.
 * @param options Further options.
 *
 * @return A sorted vector of the indices of the subsampled observations.
 */
template<typename Dim_, typename Index_, typename Value_, typename Float_>
std::vector<Index_> compute(
    Dim_ num_dims, 
    Index_ num_obs, 
    const Value_* data, 
    const knncolle::Builder<knncolle::SimpleMatrix<Dim_, Index_, Value_>, Float_>& knn_method,
    const Options& options) 
{
    auto prebuilt = knn_method.build_unique(knncolle::SimpleMatrix<Dim_, Index_, Value_>(num_dims, num_obs, data));
    return compute(*prebuilt, options);
}

}

#endif

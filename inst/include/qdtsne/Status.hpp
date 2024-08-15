/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */

#ifndef QDTSNE_STATUS_HPP
#define QDTSNE_STATUS_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <limits>

#include "SPTree.hpp"
#include "Options.hpp"
#include "utils.hpp"

/**
 * @file Status.hpp
 * @brief Status of the t-SNE algorithm.
 */

namespace qdtsne {

/**
 * @brief Current status of the t-SNE iterations.
 *
 * @tparam ndim_ Number of dimensions in the t-SNE embedding.
 * @tparam Index_ Integer type for the neighbor indices.
 * @tparam Float_ Floating-point type for the distances.
 *
 * This class holds the precomputed structures required to perform the t-SNE iterations.
 * Instances should not be constructed directly but instead created by `initialize()`.
 */
template<int ndim_, typename Index_, typename Float_> 
class Status {
public:
    /**
     * @cond
     */
    Status(NeighborList<Index_, Float_> neighbors, Options options) :
        my_neighbors(std::move(neighbors)),
        my_dY(my_neighbors.size() * ndim_), 
        my_uY(my_neighbors.size() * ndim_), 
        my_gains(my_neighbors.size() * ndim_, 1.0), 
        my_pos_f(my_neighbors.size() * ndim_), 
        my_neg_f(my_neighbors.size() * ndim_), 
        my_tree(my_neighbors.size(), options.max_depth),
        my_options(std::move(options))
    {
        if (options.num_threads > 1) {
            my_parallel_buffer.resize(my_neighbors.size());
        }
    }
    /**
     * @endcond
     */

private:
    NeighborList<Index_, Float_> my_neighbors; 
    std::vector<Float_> my_dY, my_uY, my_gains, my_pos_f, my_neg_f;

    internal::SPTree<ndim_, Float_> my_tree;
    std::vector<Float_> my_parallel_buffer; // Buffer to hold parallel-computed results prior to reduction.

    Options my_options;
    int my_iter = 0;

    typename decltype(my_tree)::LeafApproxWorkspace my_leaf_workspace;

public:
    /**
     * @return The number of iterations performed on this object so far.
     */
    int iteration() const {
        return my_iter;
    }

    /**
     * @return The maximum number of iterations. 
     * This can be modified to `run()` the algorithm for more iterations.
     */
    int max_iterations() const {
        return my_options.max_iterations;
    }

    /**
     * @return The number of observations in the dataset.
     */
    size_t num_observations() const {
        return my_neighbors.size();
    }

#ifndef NDEBUG
    /**
     * @cond
     */
    const auto& get_neighbors() const {
        return my_neighbors;
    }
    /**
     * @endcond
     */
#endif

public:
    /**
     * Run the algorithm to the specified number of iterations.
     * This can be invoked repeatedly with increasing `limit` to run the algorithm incrementally.
     *
     * @param[in, out] Y Pointer to a array containing a column-major matrix with number of rows and columns equal to `ndim_` and `num_observations()`, respectively.
     * Each row corresponds to a dimension of the embedding while each column corresponds to an observation.
     * On input, this should contain the initial location of each observation; on output, it is updated to the t-SNE location at the specified number of iterations.
     * @param limit Number of iterations to run up to.
     * The actual number of iterations performed will be the difference between `limit` and `iteration()`, i.e., `iteration()` will be equal to `limit` on completion.
     * `limit` may be greater than `max_iterations()`, to run the algorithm for more iterations than specified during construction of this `Status` object.
     */
    void run(Float_* Y, int limit) {
        Float_ multiplier = (my_iter < my_options.stop_lying_iter ? my_options.exaggeration_factor : 1);
        Float_ momentum = (my_iter < my_options.mom_switch_iter ? my_options.start_momentum : my_options.final_momentum);

        for(; my_iter < limit; ++my_iter) {
            // Stop lying about the P-values after a while, and switch momentum
            if (my_iter == my_options.stop_lying_iter) {
                multiplier = 1;
            }
            if (my_iter == my_options.mom_switch_iter) {
                momentum = my_options.final_momentum;
            }

            iterate(Y, multiplier, momentum);
        }
    }

    /**
     * Run the algorithm to the maximum number of iterations.
     * If `run()` has already been invoked with an iteration limit, this method will only perform the remaining iterations required for `iteration()` to reach `max_iter()`.
     * If `iteration()` is already greater than `max_iter()`, this method is a no-op.
     *
     * @param[in, out] Y Pointer to a array containing a column-major matrix with number of rows and columns equal to `ndim_` and `num_observations()`, respectively.
     * Each row corresponds to a dimension of the embedding while each column corresponds to an observation.
     * On input, this should contain the initial location of each observation; on output, it is updated to the t-SNE location at the specified number of iterations.
     */
    void run(Float_* Y) {
        run(Y, my_options.max_iterations);
    }

private:
    static Float_ sign(Float_ x) { 
        constexpr Float_ zero = 0;
        constexpr Float_ one = 1;
        return (x == zero ? zero : (x < zero ? -one : one));
    }

    void iterate(Float_* Y, Float_ multiplier, Float_ momentum) {
        compute_gradient(Y, multiplier);

        // Update gains
        size_t ngains = my_gains.size(); 
#ifdef _OPENMP
        #pragma omp simd
#endif
        for (size_t i = 0; i < ngains; ++i) {
            Float_& g = my_gains[i];
            constexpr Float_ lower_bound = 0.01;
            constexpr Float_ to_add = 0.2;
            constexpr Float_ to_mult = 0.8;
            g = std::max(lower_bound, sign(my_dY[i]) != sign(my_uY[i]) ? (g + to_add) : (g * to_mult));
        }

        // Perform gradient update (with momentum and gains)
#ifdef _OPENMP
        #pragma omp simd
#endif
        for (size_t i = 0; i < ngains; ++i) {
            my_uY[i] = momentum * my_uY[i] - my_options.eta * my_gains[i] * my_dY[i];
            Y[i] += my_uY[i];
        }

        // Make solution zero-mean
        size_t N = num_observations();
        for (int d = 0; d < ndim_; ++d) {
            auto start = Y + d;

            // Compute means from column-major coordinates.
            Float_ sum = 0;
            for (size_t i = 0; i < N; ++i, start += ndim_) {
                sum += *start;
            }
            sum /= N;

            start = Y + d;
            for (size_t i = 0; i < N; ++i, start += ndim_) {
                *start -= sum;
            }
        }

        return;
    }

private:
    void compute_gradient(const Float_* Y, Float_ multiplier) {
        my_tree.set(Y);
        compute_edge_forces(Y, multiplier);

        size_t N = num_observations();
        std::fill(my_neg_f.begin(), my_neg_f.end(), 0);

        Float_ sum_Q = compute_non_edge_forces();

        // Compute final t-SNE gradient
        size_t ntotal = N * static_cast<size_t>(ndim_);
        for (size_t i = 0; i < ntotal; ++i) {
            my_dY[i] = my_pos_f[i] - (my_neg_f[i] / sum_Q);
        }
    }

    void compute_edge_forces(const Float_* Y, Float_ multiplier) {
        std::fill(my_pos_f.begin(), my_pos_f.end(), 0);
        size_t N = num_observations();

#ifndef QDTSNE_CUSTOM_PARALLEL
#ifdef _OPENMP
        #pragma omp parallel num_threads(my_options.num_threads)
#endif
        {
#ifdef _OPENMP
            #pragma omp for 
#endif
            for (size_t n = 0; n < N; ++n) {
#else
        QDTSNE_CUSTOM_PARALLEL(N, [&](size_t first_, size_t last_) -> void {
            for (size_t n = first_; n < last_; ++n) {
#endif

                const auto& current = my_neighbors[n];
                size_t offset = n * static_cast<size_t>(ndim_); // cast to avoid overflow.
                const Float_* self = Y + offset;
                Float_* pos_out = my_pos_f.data() + offset;

                for (const auto& x : current) {
                    Float_ sqdist = 0; 
                    const Float_* neighbor = Y + static_cast<size_t>(x.first) * ndim_; // cast to avoid overflow.
                    for (int d = 0; d < ndim_; ++d) {
                        Float_ delta = self[d] - neighbor[d];
                        sqdist += delta * delta;
                    }

                    const Float_ mult = multiplier * x.second / (static_cast<Float_>(1) + sqdist);
                    for (int d = 0; d < ndim_; ++d) {
                        pos_out[d] += mult * (self[d] - neighbor[d]);
                    }
                }

#ifndef QDTSNE_CUSTOM_PARALLEL
            }
        }
#else
            }
        }, my_options.num_threads);
#endif

        return;
    }

    Float_ compute_non_edge_forces() {
        size_t N = num_observations();

        if (my_options.leaf_approximation) {
            my_tree.compute_non_edge_forces_for_leaves(my_options.theta, my_leaf_workspace, my_options.num_threads);
        }

#if defined(_OPENMP) || defined(QDTSNE_CUSTOM_PARALLEL)
        if (my_options.num_threads > 1) {
            // Don't use reduction methods, otherwise we get numeric imprecision
            // issues (and stochastic results) based on the order of summation.

#ifndef QDTSNE_CUSTOM_PARALLEL
#ifdef _OPENMP
            #pragma omp parallel num_threads(my_options.num_threads)
#endif
            {
#ifdef _OPENMP
                #pragma omp for
#endif
                for (size_t n = 0; n < N; ++n) {
#else
            QDTSNE_CUSTOM_PARALLEL(N, [&](size_t first_, size_t last_) -> void {
                for (size_t n = first_; n < last_; ++n) {
#endif                

                    auto neg_ptr = my_neg_f.data() + n * static_cast<size_t>(ndim_); // cast to avoid overflow.
                    if (my_options.leaf_approximation) {
                        my_parallel_buffer[n] = my_tree.compute_non_edge_forces_from_leaves(n, neg_ptr, my_leaf_workspace);
                    } else {
                        my_parallel_buffer[n] = my_tree.compute_non_edge_forces(n, my_options.theta, neg_ptr);
                    }

#ifndef QDTSNE_CUSTOM_PARALLEL
                }
            }
#else
                }
            }, my_options.num_threads);
#endif

            return std::accumulate(my_parallel_buffer.begin(), my_parallel_buffer.end(), static_cast<Float_>(0));
        }
#endif

        Float_ sum_Q = 0;
        for (size_t n = 0; n < N; ++n) {
            auto neg_ptr = my_neg_f.data() + n * static_cast<size_t>(ndim_); // cast to avoid overflow.
            if (my_options.leaf_approximation) {
                sum_Q += my_tree.compute_non_edge_forces_from_leaves(n, neg_ptr, my_leaf_workspace);
            } else {
                sum_Q += my_tree.compute_non_edge_forces(n, my_options.theta, neg_ptr);
            }
        }
        return sum_Q;
    }
};

}

#endif

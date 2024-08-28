#ifndef POWERIT_CORE_HPP
#define POWERIT_CORE_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include "aarand/aarand.hpp"

/**
 * @file core.hpp
 *
 * @brief Core data structures and calculations.
 */

#ifndef POWERIT_CUSTOM_PARALLEL
#include "subpar/subpar.hpp"
#define POWERIT_CUSTOM_PARALLEL ::subpar::parallelize
#endif

namespace powerit {

/**
 * @brief Options for `compute()`.
 */
struct Options {
    /**
     * Maximum number of iterations to perform.
     * Note that the algorithm may converge before reaching this limit.
     */
    int iterations = 500;

    /**
     * Tolerance used to determine convergence.
     * The error is defined as the L2 norm of the difference between eigenvectors at successive iterations;
     * if this drops below `tolerance`, we assume that the algorithm has converged.
     */
    double tolerance = 0.000001;

    /**
     * Number of threads to use for the matrix multiplication in `compute()`.
     * The parallelization scheme depends on the definition of the `POWERIT_CUSTOM_PARALLEL` function-like macro.
     * If undefined by the user, this macro defaults to `subpar::parallelize()` and should accept the same arguments.
     *
     * For `compute_core()`, the parallelization scheme depends on the provided `multiply()`, and this option has no effect.
     */
    int num_threads = 1;
};

/**
 * @cond
 */
template<typename Data_>
Data_ normalize(int ndim, Data_* x) {
    Data_ ss = 0;
    for (int d = 0; d < ndim; ++d) {
        ss += x[d] * x[d];
    }

    if (ss) {
        ss = std::sqrt(ss);
        for (int d = 0; d < ndim; ++d) {
            x[d] /= ss;
        }
    }
    return ss;
}
/**
 * @endcond
 */

/**
 * @brief Result of `compute()`.
 * @tparam Data_ Type of the matrix data. 
 */
template<typename Data_>
struct Result {
    /**
     * Estimate of the first eigenvalue.
     */
    Data_ value;

    /**
     * Number of iterations required to convergence.
     * Set to -1 if convergence did not occur before the maximum number of iterations specified in `Options::iterations`.
     */
    int iterations;
};

/**
 * @tparam Data_ Floating-point type for the data.
 * @tparam Engine_ Any C++11-compliant random number generator class.
 *
 * @param order Length of the array in `vector`.
 * @param[out] vector Pointer to an array of length `order`.
 * On output, this will be filled with draws from a standard normal distribution.
 * @param engine Instance of the random number generator.
 */
template<typename Data_, class Engine_>
void fill_starting_vector(size_t order, Data_* vector, Engine_& engine) {
    while (1) {
        for (size_t d = 1; d < order; d += 2) {
            auto sampled = aarand::standard_normal<Data_>(engine);
            vector[d - 1] = sampled.first;
            vector[d] = sampled.second;
        }
        if (order % 2) {
            vector[order - 1] = aarand::standard_normal<Data_>(engine).first;
        }
        if (normalize(order, vector)) {
            break;
        }
    }
}

/**
 * Perform power iterations on some diagonizable matrix to find the first eigenvalue/vector.
 * This overload generates a starting vector from an existing (P)RNG. 
 *
 * @tparam Multiply_ Function that performs a matrix dot product.
 * This is used to abstract away the representation of the actual matrix.
 * @tparam Data_ Floating-point type for the data.
 * @tparam Engine_ Any C++11-compliant random number generator class.
 *
 * @param order Order of the square matrix.
 * @param multiply Function to perform a matrix dot product with the working eigenvector.
 * This should accept `buffer`, a `std::vector<Data_>&` of length `order`;
 * and `vector`, a `const Data_*` as described below.
 * It should fill `buffer` with the product of the matrix and `vector`.
 * The return value is ignored.
 * @param[out] vector Pointer to an array of length `order`.
 * On output, this contains the estimate of the first eigenvector.
 * @param opt Further options.
 *
 * @return Result containing the first eigenvalue and other diagnostics.
 */
template<class Multiply_, typename Data_>
Result<Data_> compute_core(size_t order, Multiply_ multiply, Data_* vector, const Options& opt) {
    Result<Data_> stats;
    auto& l2 = stats.value;
    stats.iterations = -1;
    std::vector<Data_> buffer(order);

    for (int i = 0; i < opt.iterations; ++i) {
        multiply(buffer, vector);
        l2 = normalize(order, buffer.data());

        // Assuming convergence if the vector did not change much from the last iteration.
        Data_ err = 0;
        for (size_t d = 0; d < order; ++d) {
            Data_ diff = buffer[d] - vector[d];
            err += diff * diff;
        }
        if (std::sqrt(err) < opt.tolerance) {
            stats.iterations = i + 1;
            break;
        }

        std::copy(buffer.begin(), buffer.end(), vector);
    }

    return stats;
} 

}

#endif

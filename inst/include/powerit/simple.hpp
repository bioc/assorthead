#ifndef POWERIT_SIMPLE_HPP
#define POWERIT_SIMPLE_HPP

#include "core.hpp"
#include <numeric>

/**
 * @file simple.hpp
 *
 * @brief Power iterations for a simple array.
 */

namespace powerit {

/**
 * Perform power iterations on a diagonizable matrix to find the first eigenvalue/vector.
 * This overload generates a starting vector from an existing (P)RNG. 
 *
 * @tparam Data_ Floating-point type for the data.
 * @tparam Engine_ Any C++11-compliant random number generator class.
 *
 * @param order Order of the square matrix.
 * @param[in] matrix Pointer to an array containing an `order`-by-`order` diagonalizable matrix.
 * @param row_major Whether `matrix` is row-major.
 * @param[out] vector Pointer to an array of length `order`.
 * On output, this contains the estimate of the first eigenvector.
 * @param engine Instance of the random number generator.
 * @param opt Further options.
 *
 * @return Result containing the first eigenvalue and other diagnostics.
 */
template<typename Data_, class Engine_>
Result<Data_> compute(size_t order, const Data_* matrix, bool row_major, Data_* vector, Engine_& engine, const Options& opt) {
    fill_starting_vector(order, vector, engine);
    return compute(order, matrix, row_major, vector, opt);
}

/**
 * Perform power iterations on an array containing a diagonizable matrix. 
 * This overload assumes that a random starting vector has already been generated.
 *
 * @tparam Data_ Floating-point type for the data.
 *
 * @param order Order of the square matrix.
 * @param[in] matrix Pointer to an array containing an `order`-by-`order` diagonalizable matrix.
 * @param row_major Whether `matrix` is row-major.
 * @param[in,out] vector Pointer to an array of length `order`.
 * On input, this should contain a random starting vector.
 * On output, this contains the estimate of the first eigenvector.
 * @param opt Further options.
 *
 * @return Result containing the first eigenvalue and other diagnostics.
 */
template<typename Data_>
Result<Data_> compute(size_t order, const Data_* matrix, bool row_major, Data_* vector, const Options& opt) {
    if (row_major) {
        return compute_core(order, [&](std::vector<Data_>& buffer, const Data_* vec) {
            POWERIT_CUSTOM_PARALLEL(opt.num_threads, order, [&](int, size_t start, size_t length) {
                for (size_t j = start, end = start + length; j < end; ++j) {
                    // Note that j and order are already both 'size_t', so no need to cast to avoid overflow.
                    buffer[j] = std::inner_product(vec, vec + order, matrix + j * order, static_cast<Data_>(0.0));
                }
            });
        }, vector, opt);

    } else if (opt.num_threads == 1) { 
        // Dedicated path to avoid allocating a per-thread temporary.
        return compute_core(order, [&](std::vector<Data_>& buffer, const Data_* vec) {
            std::fill(buffer.begin(), buffer.end(), 0);
            auto matcopy = matrix;
            for (size_t j = 0; j < order; ++j) {
                Data_ mult = vec[j];
                for (size_t k = 0; k < order; ++k, ++matcopy) {
                    buffer[k] += mult * (*matcopy);
                }
            }
        }, vector, opt);

    } else {
        return compute_core(order, [&](std::vector<Data_>& buffer, const Data_* vec) {
            POWERIT_CUSTOM_PARALLEL(opt.num_threads, order, [&](int, size_t start, size_t length) {
                std::vector<Data_> tmp(length);
                size_t offset = start; // already size_t's, no need to cast.
                for (size_t j = 0; j < order; ++j, offset += order) {
                    auto mult = vec[j];
                    auto matcopy = matrix + offset;
                    for (size_t k = 0; k < length; ++k, ++matcopy) {
                        tmp[k] += mult * (*matcopy);
                    }
                }
                std::copy(tmp.begin(), tmp.end(), buffer.begin() + start);
            });
        }, vector, opt);
    } 
}

}

#endif

#ifndef TATAMI_STATS_GROUPED_VARIANCES_HPP
#define TATAMI_STATS_GROUPED_VARIANCES_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include "variances.hpp"
#include <vector>
#include <algorithm>

/**
 * @file grouped_variances.hpp
 *
 * @brief Compute group-wise variances from a `tatami::Matrix`.
 */

namespace tatami_stats {

/**
 * @brief Functions for computing dimension-wise grouped variances.
 * @namespace tatami_stats::grouped_variances
 */
namespace grouped_variances {

/**
 * @brief Grouped summation options.
 */
struct Options {
    /**
     * Whether to check for NaNs in the input, and skip them.
     * If false, NaNs are assumed to be absent, and the behavior of the summation in the presence of NaNs is undefined.
     */
    bool skip_nan = false;

    /**
     * Number of threads to use when computing variances across a `tatami::Matrix`.
     * See `tatami::parallelize()` for more details on the parallelization mechanism.
     */
    int num_threads = 1;
};

/**
 * @cond
 */
namespace internal {

template<typename Index_, typename Output_>
void finish_means(size_t num_groups, const Index_* group_size, Output_* output_means) {
    for (size_t b = 0; b < num_groups; ++b) {
        if (group_size[b]) {
            output_means[b] /= group_size[b];
        } else {
            output_means[b] = std::numeric_limits<Output_>::quiet_NaN();
        }
    }
}

template<typename Index_, typename Output_>
void finish_variances(size_t num_groups, const Index_* group_size, Output_* output_variances) {
    for (size_t b = 0; b < num_groups; ++b) {
        if (group_size[b] > 1) {
            output_variances[b] /= group_size[b] - 1;
        } else {
            output_variances[b] = std::numeric_limits<Output_>::quiet_NaN();
        }
    }
}

}
/**
 * @endcond
 */

/**
 * Compute the mean and variance from a dense objective vector.
 * This uses the standard two-pass algorithm with naive accumulation of the sum of squared differences;
 * thus, it is best used with a sufficiently high-precision `Output_` like `double`.
 *
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Integer type of the matrix indices.
 * @tparam Group_ Integer type of the group assignments.
 * @tparam Output_ Type of the output data.
 *
 * @param[in] ptr Pointer to an array of values of length `num`.
 * @param num Length of the objective vector, i.e., length of the array at `ptr`.
 * @param[in] group Pointer to an array of length `num`, containing the group assignment for each entry of `ptr`.
 * Entries of `group` should lie in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param num_groups Number of groups, i.e., \f$N\f$.
 * @param[in] group_size Pointer to an array of length `num_groups`, containing the size of each group.
 * This can be obtained by calling `tabulate_groups()` on `group`.
 * @param[out] output_means Pointer to an array of length `num_groups`.
 * This is filled with the per-group mean on output.
 * Values may be NaN if there are not enough (non-NaN) values in a group.
 * @param[out] output_variances Pointer to an array of length `num_groups`.
 * This is filled with the per-group variances on output.
 * Values may be NaN if there are not enough (non-NaN) values in a group.
 * @param skip_nan See `Options::skip_nan`.
 * @param[out] valid_group_size Pointer to an array of length `num_groups`.
 * This is used to store the number of non-NaN entries.
 * Only used if `skip_nan = true`.
 */
template<typename Value_, typename Index_, typename Group_, typename Output_>
void direct(
    const Value_* ptr, 
    Index_ num, 
    const Group_* group, 
    size_t num_groups, 
    const Index_* group_size, 
    Output_* output_means, 
    Output_* output_variances, 
    bool skip_nan, 
    Index_* valid_group_size)
{
    std::fill_n(output_means, num_groups, 0);
    std::fill_n(output_variances, num_groups, 0);

    ::tatami_stats::internal::nanable_ifelse<Value_>(
        skip_nan,
        [&]() {
            std::fill_n(valid_group_size, num_groups, 0);

            for (Index_ j = 0; j < num; ++j) {
                auto x = ptr[j];
                if (!std::isnan(x)) {
                    auto b = group[j];
                    output_means[b] += x;
                    ++valid_group_size[b];
                }
            }
            internal::finish_means(num_groups, valid_group_size, output_means);

            for (Index_ j = 0; j < num; ++j) {
                auto x = ptr[j];
                if (!std::isnan(x)) {
                    auto b = group[j];
                    auto delta = x - output_means[b];
                    output_variances[b] += delta * delta;
                }
            }
            internal::finish_variances(num_groups, valid_group_size, output_variances);
        },
        [&]() {
            for (Index_ j = 0; j < num; ++j) {
                output_means[group[j]] += ptr[j];
            }
            internal::finish_means(num_groups, group_size, output_means);

            for (Index_ j = 0; j < num; ++j) {
                auto b = group[j];
                auto delta = ptr[j] - output_means[b];
                output_variances[b] += delta * delta;
            }
            internal::finish_variances(num_groups, group_size, output_variances);
        }
    );
}

/**
 * Compute the mean and variance from a sparse objective vector.
 * This uses the standard two-pass algorithm with naive accumulation of the sum of squared differences;
 * thus, it is best used with a sufficiently high-precision `Output_` like `double`.
 *
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Integer type of the matrix indices.
 * @tparam Group_ Integer type of the group assignments.
 * @tparam Output_ Type of the output data.
 *
 * @param[in] value Pointer to an array of length `num_nonzero`, containing the values of the structural non-zeros.
 * @param[in] index Pointer to an array of length `num_nonzero`, containing the indices of the structural non-zeros.
 * All indices should be non-negative and less than the length of the objective vector.
 * @param num_nonzero Number of structural non-zeros. 
 * @param[in] group Pointer to an array of length equal to the length of the objective vector, containing the group assignment for each vector element.
 * Entries of `group` should lie in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param num_groups Number of groups, i.e., \f$N\f$.
 * @param[in] group_size Pointer to an array of length `num_groups`, containing the size of each group.
 * This can be obtained by calling `tabulate_groups()` on `group`.
 * @param[out] output_means Pointer to an array of length `num_groups`.
 * This is filled with the per-group mean on output.
 * Values may be NaN if there are not enough (non-NaN) values in a group.
 * @param[out] output_variances Pointer to an array of length `num_groups`.
 * This is filled with the per-group variances on output.
 * Values may be NaN if there are not enough (non-NaN) values in a group.
 * @param[out] output_nonzero Pointer to an array of length `num_groups`.
 * On output, this is filled with the number of structural non-zeros in each group.
 * @param skip_nan See `Options::skip_nan`.
 * @param[out] valid_group_size Pointer to an array of length `num_groups`.
 * This is used to store the number of non-NaN entries.
 * Only used if `skip_nan = true`.
 */
template<typename Value_, typename Index_, typename Group_, typename Output_>
void direct(
    const Value_* value, 
    const Index_* index, 
    Index_ num_nonzero, 
    const Group_* group, 
    size_t num_groups, 
    const Index_* group_size, 
    Output_* output_means, 
    Output_* output_variances, 
    Index_* output_nonzero,
    bool skip_nan, 
    Index_* valid_group_size)
{
    std::fill_n(output_means, num_groups, 0);
    std::fill_n(output_nonzero, num_groups, 0);
    std::fill_n(output_variances, num_groups, 0);

    ::tatami_stats::internal::nanable_ifelse<Value_>(
        skip_nan,
        [&]() {
            std::copy_n(group_size, num_groups, valid_group_size);

            for (Index_ j = 0; j < num_nonzero; ++j) {
                auto x = value[j];
                auto b = group[index[j]];
                if (!std::isnan(x)) {
                    output_means[b] += x;
                    ++(output_nonzero[b]);
                } else {
                    --(valid_group_size[b]);
                }
            }
            internal::finish_means(num_groups, valid_group_size, output_means);

            for (Index_ j = 0; j < num_nonzero; ++j) {
                auto x = value[j];
                if (!std::isnan(x)) {
                    auto b = group[index[j]];
                    auto delta = x - output_means[b];
                    output_variances[b] += delta * delta;
                }
            }
            for (size_t b = 0; b < num_groups; ++b) {
                output_variances[b] += output_means[b] * output_means[b] * (valid_group_size[b] - output_nonzero[b]);
            }
            internal::finish_variances(num_groups, valid_group_size, output_variances);
        },
        [&]() {
            for (Index_ j = 0; j < num_nonzero; ++j) {
                auto b = group[index[j]];
                output_means[b] += value[j];
                ++output_nonzero[b];
            }
            internal::finish_means(num_groups, group_size, output_means);

            for (Index_ j = 0; j < num_nonzero; ++j) {
                auto b = group[index[j]];
                auto delta = value[j] - output_means[b];
                output_variances[b] += delta * delta;
            }
            for (size_t b = 0; b < num_groups; ++b) {
                output_variances[b] += output_means[b] * output_means[b] * (group_size[b] - output_nonzero[b]);
            }
            internal::finish_variances(num_groups, group_size, output_variances);
        }
    );
}

/**
 * Compute per-group variances for each element of a chosen dimension of a `tatami::Matrix`.
 *
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 * @tparam Output_ Type of the output value.
 * This should be floating-point to store potential averages.
 *
 * @param row Whether to compute variances for the rows.
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns (if `row = true`) or rows (otherwise).
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param num_groups Number of groups, i.e., \f$N\f$.
 * This can be determined by calling `tatami_stats::total_groups()` on `group`.
 * @param[in] group_size Pointer to an array of length `num_groups`, containing the size of each group.
 * @param[out] output Pointer to an array of pointers of length equal to the number of groups.
 * Each inner pointer should reference an array of length equal to the number of rows (if `row = true`) or columns (otherwise).
 * On output, this will contain the row/column variances for each group (indexed according to the assignment in `group`).
 * @param sopt Summation options.
 */
template<typename Value_, typename Index_, typename Group_, typename Output_>
void apply(bool row, const tatami::Matrix<Value_, Index_>* p, const Group_* group, size_t num_groups, const Index_* group_size, Output_** output, const Options& sopt) {
    Index_ dim = (row ? p->nrow() : p->ncol());
    Index_ otherdim = (row ? p->ncol() : p->nrow());

    if (p->sparse()) {
        if (p->prefer_rows() == row) {
            tatami::parallelize([&](int, Index_ start, Index_ len) -> void {
                auto ext = tatami::consecutive_extractor<true>(p, row, start, len);
                std::vector<Value_> xbuffer(otherdim);
                std::vector<Index_> ibuffer(otherdim);

                std::vector<Output_> tmp_means(num_groups);
                std::vector<Output_> output_variances(num_groups);
                std::vector<Index_> tmp_nonzero(num_groups);
                std::vector<Index_> valid_group_size(sopt.skip_nan ? num_groups : 0);

                for (Index_ i = 0; i < len; ++i) {
                    auto range = ext->fetch(xbuffer.data(), ibuffer.data());
                    direct(
                        range.value, 
                        range.index, 
                        range.number, 
                        group, 
                        num_groups, 
                        group_size, 
                        tmp_means.data(), 
                        output_variances.data(), 
                        tmp_nonzero.data(), 
                        sopt.skip_nan, 
                        valid_group_size.data()
                    );

                    for (size_t g = 0; g < num_groups; ++g) {
                        output[g][i + start] = output_variances[g];
                    }
                }
            }, dim, sopt.num_threads);

        } else {
            // Order within each observed vector doesn't affect numerical
            // precision of the outcome, as addition order for each objective
            // vector is already well-defined for a running calculation.
            tatami::Options opt;
            opt.sparse_ordered_index = false; 

            tatami::parallelize([&](int thread, Index_ start, Index_ len) -> void {
                std::vector<variances::RunningSparse<Output_, Value_, Index_> > runners;
                runners.reserve(num_groups);
                std::vector<LocalOutputBuffer<Output_> > local_var_output;
                local_var_output.reserve(num_groups);
                std::vector<std::vector<Output_> > local_mean_output;
                local_mean_output.reserve(num_groups);

                for (size_t g = 0; g < num_groups; ++g) {
                    local_var_output.emplace_back(thread, start, len, output[g]);
                    local_mean_output.emplace_back(len);
                    runners.emplace_back(len, local_mean_output.back().data(), local_var_output.back().data(), sopt.skip_nan, start);
                }

                auto ext = tatami::consecutive_extractor<true>(p, !row, static_cast<Index_>(0), otherdim, start, len, opt);
                std::vector<Value_> xbuffer(len);
                std::vector<Index_> ibuffer(len);

                for (int i = 0; i < otherdim; ++i) {
                    auto range = ext->fetch(xbuffer.data(), ibuffer.data());
                    runners[group[i]].add(range.value, range.index, range.number);
                }

                for (size_t g = 0; g < num_groups; ++g) {
                    runners[g].finish();
                    local_var_output[g].transfer();
                }
            }, dim, sopt.num_threads);
        }

    } else {
        if (p->prefer_rows() == row) {
            tatami::parallelize([&](int, Index_ start, Index_ len) -> void {
                auto ext = tatami::consecutive_extractor<false>(p, row, start, len);
                std::vector<Value_> xbuffer(otherdim);
                std::vector<Output_> tmp_means(num_groups);
                std::vector<Output_> output_variances(num_groups);
                std::vector<Index_> valid_group_size(sopt.skip_nan ? num_groups : 0);

                for (Index_ i = 0; i < len; ++i) {
                    auto ptr = ext->fetch(xbuffer.data());
                    direct(
                        ptr, 
                        otherdim, 
                        group, 
                        num_groups, 
                        group_size, 
                        tmp_means.data(), 
                        output_variances.data(), 
                        sopt.skip_nan, 
                        valid_group_size.data()
                    );

                    for (size_t g = 0; g < num_groups; ++g) {
                        output[g][i + start] = output_variances[g];
                    }
                }
            }, dim, sopt.num_threads);

        } else {
            tatami::parallelize([&](int thread, Index_ start, Index_ len) -> void {
                std::vector<variances::RunningDense<Output_, Value_, Index_> > runners;
                runners.reserve(num_groups);
                std::vector<LocalOutputBuffer<Output_> > local_var_output;
                local_var_output.reserve(num_groups);
                std::vector<std::vector<Output_> > local_mean_output;
                local_mean_output.reserve(num_groups);

                for (size_t g = 0; g < num_groups; ++g) {
                    local_var_output.emplace_back(thread, start, len, output[g]);
                    local_mean_output.emplace_back(len);
                    runners.emplace_back(len, local_mean_output.back().data(), local_var_output.back().data(), sopt.skip_nan);
                }

                std::vector<Value_> xbuffer(len);
                auto ext = tatami::consecutive_extractor<false>(p, !row, static_cast<Index_>(0), otherdim, start, len);

                for (Index_ i = 0; i < otherdim; ++i) {
                    auto ptr = ext->fetch(xbuffer.data());
                    runners[group[i]].add(ptr);
                }

                for (size_t g = 0; g < num_groups; ++g) {
                    runners[g].finish();
                    local_var_output[g].transfer();
                }
            }, dim, sopt.num_threads);
        }
    }
}

/**
 * Wrapper around `apply()` for row-wise grouped variances.
 *
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each row.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param sopt Summation options.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of rows, containing the row-wise variances for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_row(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const Options& sopt) {
    size_t mydim = p->nrow();
    auto group_size = tabulate_groups(group, p->ncol());
    size_t ngroup = group_size.size();

    std::vector<std::vector<Output_> > output(ngroup);
    std::vector<Output_*> ptrs;
    ptrs.reserve(output.size());
    for (auto& o : output) {
        o.resize(mydim);
        ptrs.push_back(o.data());
    }

    apply(true, p, group, ngroup, group_size.data(), ptrs.data(), sopt);
    return output;
}

/**
 * Overload with default options.
 *
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of rows, containing the row-wise variances for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_row(const tatami::Matrix<Value_, Index_>* p, const Group_* group) {
    return by_row(p, group, Options());
}

/**
 * Wrapper around `apply()` for column-wise grouped variances.
 *
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the column/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of rows.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 * @param sopt Summation options.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of columns, containing the column-wise variances for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_column(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const Options& sopt) {
    size_t mydim = p->ncol();
    auto group_size = tabulate_groups(group, p->nrow());
    size_t ngroup = group_size.size();

    std::vector<std::vector<Output_> > output(ngroup);
    std::vector<Output_*> ptrs;
    ptrs.reserve(output.size());
    for (auto& o : output) {
        o.resize(mydim);
        ptrs.push_back(o.data());
    }

    apply(false, p, group, ngroup, group_size.data(), ptrs.data(), sopt);
    return output;
}

/**
 * Overload with default options.
 *
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the column/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of rows.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in \f$[0, N)\f$ where \f$N\f$ is the number of unique groups.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of columns, containing the column-wise variances for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_column(const tatami::Matrix<Value_, Index_>* p, const Group_* group) {
    return by_column(p, group, Options());
}

}

}

#endif

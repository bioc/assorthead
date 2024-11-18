#ifndef SCRAN_VARIANCES_CHOOSE_HIGHLY_VARIABLE_GENES_HPP
#define SCRAN_VARIANCES_CHOOSE_HIGHLY_VARIABLE_GENES_HPP

#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdint>

/**
 * @file choose_highly_variable_genes.hpp
 * @brief Choose highly variable genes for downstream analyses.
 */

namespace scran_variances {

/**
 * @brief Options for `choose_highly_variable_genes()`.
 */
struct ChooseHighlyVariableGenesOptions {
    /**
     * Number of top genes to choose.
     * Note that the actual number of chosen genes may be: 
     *
     * - smaller than `top`, if the latter is greater than the total number of genes in the dataset. 
     * - smaller than `top`, if `ChooseHighlyVariableGenesOptions::use_bound = true` and `top` is greater than the total number of genes in the dataset with statistics greater than `ChooseHighlyVariableGenesOptions::bound`.
     * - larger than `top`, if `ChooseHighlyVariableGenesOptions::keep_ties = true` and there are multiple ties at the `top`-th chosen gene.
     */
    size_t top = 4000;

    /**
     * Whether larger statistics correspond to higher variances.
     */
    bool larger = true;

    /**
     * Whether to consider an absolute bound on the statistic when choosing HVGs.
     * The value of the bound is determined by `ChooseHighlyVariableGenesOptions::bound`.
     */
    bool use_bound = false;

    /**
     * A lower bound for the statistic, at or below which a gene will not be considered as highly variable even if it is among the top `top` genes.
     * If `ChooseHighlyVariableGenesOptions::larger = false`, this is an upper bound instead.
     * Only used if `ChooseHighlyVariableGenesOptions::use_bound = true`.
     */
    double bound = 0;

    /**
     * Whether to keep all genes with statistics that are tied with the `ChooseHighlyVariableGenesOptions::top`-th gene.
     * If `false`, ties are arbitrarily broken but the number of retained genes will not be greater than `ChooseHighlyVariableGenesOptions::top`.
     */
    bool keep_ties = true;
};

/**
 * @cond
 */
namespace internal {

template<typename Index_, typename Stat_, class Cmp_>
std::vector<Index_> create_semisorted_indices(size_t n, const Stat_* statistic, Cmp_ cmp, size_t top) {
    std::vector<Index_> collected(n);
    std::iota(collected.begin(), collected.end(), static_cast<Index_>(0));
    auto cBegin = collected.begin(), cMid = cBegin + top - 1, cEnd = collected.end();
    std::nth_element(cBegin, cMid, cEnd, [&](Index_ l, Index_ r) -> bool { 
        auto L = statistic[l], R = statistic[r];
        if (L == R) {
            return l < r; // always favor the earlier index for a stable sort, even if options.larger = false.
        } else {
            return cmp(L, R);
        }
    });
    return collected;
}

template<typename Stat_, class Output_, class Cmp_, class CmpEqual_>
void choose_highly_variable_genes(size_t n, const Stat_* statistic, Output_* output, Cmp_ cmp, CmpEqual_ cmpeq, const ChooseHighlyVariableGenesOptions& options) {
    if (options.top == 0) {
        std::fill_n(output, n, false);
        return;
    }

    Stat_ bound = options.bound;
    if (options.top >= n) {
        if (options.use_bound) {
            for (size_t i = 0; i < n; ++i) {
                output[i] = cmp(statistic[i], bound);
            }
        } else {
            std::fill_n(output, n, true);
        }
        return;
    }

    auto collected = create_semisorted_indices<size_t>(n, statistic, cmp, options.top);
    Stat_ threshold = statistic[collected[options.top - 1]];

    if (options.keep_ties) {
        if (options.use_bound && !cmp(threshold, bound)) {
            for (size_t i = 0; i < n; ++i) {
                output[i] = cmp(statistic[i], bound);
            }
        } else {
            for (size_t i = 0; i < n; ++i) {
                output[i] = cmpeq(statistic[i], threshold);
            }
        }
        return;
    }

    std::fill_n(output, n, false);
    size_t counter = options.top;
    if (options.use_bound && !cmp(threshold, bound)) {
        --counter;
        while (counter > 0) {
            --counter;
            if (cmp(statistic[collected[counter]], bound)) {
                ++counter;
                break;
            }
        }
    }

    for (size_t i = 0; i < counter; ++i) {
        output[collected[i]] = true;
    }
}

template<typename Index_, typename Stat_, class Cmp_, class CmpEqual_>
std::vector<Index_> choose_highly_variable_genes_index(size_t n, const Stat_* statistic, Cmp_ cmp, CmpEqual_ cmpeq, const ChooseHighlyVariableGenesOptions& options) {
    std::vector<Index_> output;
    if (options.top == 0) {
        return output;
    }

    Stat_ bound = options.bound;
    if (options.top >= n) {
        if (options.use_bound) {
            for (size_t i = 0; i < n; ++i) {
                if (options.use_bound && cmp(statistic[i], bound)) {
                    output.push_back(i);
                }
            }
        } else {
            output.resize(n);
            std::iota(output.begin(), output.end(), static_cast<Index_>(0));
        }
        return output;
    }

    output = create_semisorted_indices<Index_>(n, statistic, cmp, options.top);
    Stat_ threshold = statistic[output[options.top - 1]];

    if (options.keep_ties) {
        output.clear();
        if (options.use_bound && !cmp(threshold, bound)) {
            for (size_t i = 0; i < n; ++i) {
                if (cmp(statistic[i], bound)) {
                    output.push_back(i);
                }
            }
        } else {
            for (size_t i = 0; i < n; ++i) {
                if (cmpeq(statistic[i], threshold)) {
                    output.push_back(i);
                }
            }
        }
        return output;
    }

    size_t counter = options.top;
    if (options.use_bound && !cmp(threshold, bound)) {
        --counter;
        while (counter > 0) {
            --counter;
            if (cmp(statistic[output[counter]], bound)) {
                ++counter;
                break;
            }
        }
    }

    output.resize(counter);
    std::sort(output.begin(), output.end());
    return output;
}

}
/**
 * @endcond
 */

/**
 * @tparam Stat_ Type of the variance statistic.
 * @tparam Bool_ Type to be used as a boolean.
 *
 * @param n Number of genes.
 * @param[in] statistic Pointer to an array of length `n` containing the per-gene variance statistics.
 * @param[out] output Pointer to an array of length `n`. 
 * On output, this is filled with `true` if the gene is to be retained and `false` otherwise.
 * @param options Further options.
 */
template<typename Stat_, typename Bool_>
void choose_highly_variable_genes(size_t n, const Stat_* statistic, Bool_* output, const ChooseHighlyVariableGenesOptions& options) {
    if (options.larger) {
        internal::choose_highly_variable_genes(
            n, 
            statistic, 
            output, 
            [](Stat_ l, Stat_ r) -> bool { return l > r; },
            [](Stat_ l, Stat_ r) -> bool { return l >= r; },
            options
        );
    } else {
        internal::choose_highly_variable_genes(
            n, 
            statistic, 
            output, 
            [](Stat_ l, Stat_ r) -> bool { return l < r; },
            [](Stat_ l, Stat_ r) -> bool { return l <= r; },
            options
        );
    }
}

/**
 * @tparam Stat_ Type of the variance statistic.
 * @tparam Bool_ Type to be used as a boolean.
 *
 * @param n Number of genes.
 * @param[in] statistic Pointer to an array of length `n` containing the per-gene variance statistics.
 * @param options Further options.
 *
 * @return A vector of booleans of length `n`, indicating whether each gene is to be retained.
 */
template<typename Bool_ = uint8_t, typename Stat_>
std::vector<Bool_> choose_highly_variable_genes(size_t n, const Stat_* statistic, const ChooseHighlyVariableGenesOptions& options) {
    std::vector<Bool_> output(n);
    choose_highly_variable_genes(n, statistic, output.data(), options);
    return output;
}

/**
 * @tparam Index_ Type of the indices.
 * @tparam Stat_ Type of the variance statistic.
 *
 * @param n Number of genes.
 * @param[in] statistic Pointer to an array of length `n` containing the per-gene variance statistics.
 * @param options Further options.
 *
 * @return Vector of sorted and unique indices for the chosen genes.
 * All indices are guaranteed to be non-negative and less than `n`.
 */
template<typename Index_, typename Stat_>
std::vector<Index_> choose_highly_variable_genes_index(Index_ n, const Stat_* statistic, const ChooseHighlyVariableGenesOptions& options) {
    if (options.larger) {
        return internal::choose_highly_variable_genes_index<Index_>(
            n, 
            statistic, 
            [](Stat_ l, Stat_ r) -> bool { return l > r; },
            [](Stat_ l, Stat_ r) -> bool { return l >= r; },
            options
        );
    } else {
        return internal::choose_highly_variable_genes_index<Index_>(
            n, 
            statistic, 
            [](Stat_ l, Stat_ r) -> bool { return l < r; },
            [](Stat_ l, Stat_ r) -> bool { return l <= r; },
            options
        );
    }
}

}

#endif

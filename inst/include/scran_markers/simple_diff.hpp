#ifndef SCRAN_MARKERS_SIMPLE_DIFF_HPP
#define SCRAN_MARKERS_SIMPLE_DIFF_HPP

#include <limits>

#include "PrecomputedPairwiseWeights.hpp"

namespace scran_markers {

namespace internal {

// 'values' is expected to be an 'ngroups * nblocks' array where groups are the
// faster-changing dimension and the blocks are slower.
template<typename Stat_, typename Weight_>
Stat_ compute_pairwise_simple_diff(size_t g1, size_t g2, const Stat_* values, size_t ngroups, size_t nblocks, const PrecomputedPairwiseWeights<Weight_>& preweights) {
    auto winfo = preweights.get(g1, g2);
    auto total_weight = winfo.second;
    if (total_weight == 0) {
        return std::numeric_limits<Stat_>::quiet_NaN();
    }

    Stat_ output = 0;
    size_t offset1 = g1, offset2 = g2; 
    for (size_t b = 0; b < nblocks; ++b, offset1 += ngroups, offset2 += ngroups) {
        auto weight = winfo.first[b];
        if (weight) {
            auto left = values[offset1 /* == b * ngroups + g1 */];
            auto right = values[offset2 /* == b * ngroups + g2 */]; 
            output += (left - right) * weight;
        }
    }

    return output / total_weight;
}

template<typename Stat_, typename Weight_>
void compute_pairwise_simple_diff(const Stat_* values, size_t ngroups, size_t nblocks, const PrecomputedPairwiseWeights<Weight_>& preweights, Stat_* output) {
    size_t offset1_raw = ngroups, offset2_raw = 1;
    for (size_t g1 = 1; g1 < ngroups; ++g1, offset1_raw += ngroups, ++offset2_raw) {
        auto offset1 = offset1_raw, offset2 = offset2_raw;
        for (size_t g2 = 0; g2 < g1; ++g2, ++offset1, offset2 += ngroups) {
            auto d = compute_pairwise_simple_diff(g1, g2, values, ngroups, nblocks, preweights);
            output[offset1 /* == g1 * ngroups + g2 */] = d;
            output[offset2 /* == g2 * ngroups + g1 */] = -d;
        }
    }
}

}

}

#endif

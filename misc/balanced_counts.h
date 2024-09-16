#ifndef GET_BALANCED_COUNTS_H_
#define GET_BALANCED_COUNTS_H_

#include <vector> // std::vector
#include <concepts> // std::integral
#include <algorithm> // std::fill
#include "mytypeinfo.h"

template <index_type Index>
void get_balanced_counts(std::vector<Index>& counts, size_t totsize)
{
    Index blocks = counts.size();
    std::fill(counts.begin(), counts.end(), totsize/blocks);

    counts.back() = totsize - (blocks-1)*(totsize/blocks);

    Index diff = counts.back() - counts.front();

    for (Index i = 0; i < diff; ++i)
    {
        counts[blocks-1-i]++;
        counts[blocks-1]--;
    }
}

#endif

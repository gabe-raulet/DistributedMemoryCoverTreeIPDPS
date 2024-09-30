#ifndef INSERT_TREE_H_
#define INSERT_TREE_H_

#include "mytypeinfo.h"
#include <numeric>
#include <vector>
#include <string>
#include <sstream>
#include "json.hpp"

using json = nlohmann::json;

template <class Item_, index_type Index_>
struct InsertTree
{
    using Item = Item_;
    using Index = Index_;

    using IndexVector = std::vector<Index>;

    InsertTree() : nlevels(0) {}

    Index add_vertex(Item item, Index parent);

    Index get_children(Index parent, IndexVector& ids) const;
    Index num_children(Index parent) const { return child_displs[parent+1]-child_displs[parent]; }
    bool is_leaf(Index id) const { return num_children(id)==0; }

    Item operator[](Index id) const { return vertices[id]; }
    Item& operator[](Index id) { return vertices[id]; }

    Index num_levels() const { return nlevels; }
    Index num_vertices() const { return vertices.size(); }

    void clear();

    using ItemVector = std::vector<Item>;
    using IndexVectorVector = std::vector<IndexVector>;

    template <class Itemizer>
    void get_json_repr(json& json_repr, Itemizer itemizer) const;

    template <class Itemizer>
    std::string get_json_repr(Itemizer itemizer) const
    {
        json json_repr;
        get_json_repr(json_repr, itemizer);
        std::stringstream ss;
        ss << std::setw(4) << json_repr << std::endl;
        return ss.str();
    }

    ItemVector vertices;
    IndexVector levels;
    IndexVector parents;
    IndexVector child_displs;
    IndexVector child_counts;
    IndexVector child_vals;
    Index nlevels;

    void fill_child_vals()
    {
        Index n = num_vertices();
        child_displs.resize(n+1, 0);

        for (Index i = 1; i <= n; ++i)
        {
            child_displs[i] = child_displs[i-1] + child_counts[i-1];
        }

        IndexVector child_ptrs = child_displs;

        child_vals.resize(child_displs.back());

        for (Index i = 1; i < n; ++i)
        {
            Index parent = parents[i];
            Index loc = child_ptrs[parent]++;
            child_vals[loc] = i;
        }
    }
};

#include "itree.hpp"

#endif

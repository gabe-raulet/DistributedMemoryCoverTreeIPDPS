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
    Index num_children(Index parent) const { return children[parent].size(); }
    bool is_leaf(Index id) const { return children[id].empty(); }

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
    IndexVectorVector children;
    Index nlevels;
};

#include "itree.hpp"

#endif

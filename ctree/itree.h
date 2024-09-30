#ifndef INSERT_TREE_H_
#define INSERT_TREE_H_

#include "mytypeinfo.h"
#include <numeric>
#include <vector>
#include <string>
#include <sstream>
#include <omp.h>
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

    template <class Itemizer, class NewItem>
    void itemize_new_tree(InsertTree<NewItem, Index>& new_tree, Itemizer itemizer, bool threaded=false) const
    {
        new_tree.clear();
        new_tree.levels = levels;
        new_tree.parents = parents;
        new_tree.children = children;
        new_tree.nlevels = nlevels;

        Index n = num_vertices();
        new_tree.vertices.resize(n);

        #pragma omp parallel for if (threaded)
        for (Index i = 0; i < n; ++i)
        {
            new_tree.vertices[i] = itemizer(vertices[i]);
        }
    }

    ItemVector vertices;
    IndexVector levels;
    IndexVector parents;
    IndexVectorVector children;
    Index nlevels;
};

#include "itree.hpp"

#endif

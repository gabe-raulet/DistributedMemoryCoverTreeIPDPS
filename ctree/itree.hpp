template <class Item_, index_type Index_>
typename InsertTree<Item_, Index_>::Index
InsertTree<Item_, Index_>::add_vertex(Item item, Index parent)
{
    Index level;
    Index vertex = vertices.size();

    vertices.push_back(item);
    parents.push_back(parent);
    children.emplace_back();

    if (parent >= 0)
    {
        level = levels[parent] + 1;
        children[parent].push_back(vertex);
    }
    else level = 0;

    nlevels = std::max(level+1, nlevels);
    levels.push_back(level);

    return vertex;
}

template <class Item_, index_type Index_>
typename InsertTree<Item_, Index_>::Index
InsertTree<Item_, Index_>::get_children(Index parent, IndexVector& ids) const
{
    ids = children[parent];
    return ids.size();
}

template <class Item_, index_type Index_>
void InsertTree<Item_, Index_>::clear()
{
    parents.clear();
    vertices.clear();
    levels.clear();
    children.clear();
    nlevels = 0;
}

template <class Item_, index_type Index_>
template <class Itemizer>
void InsertTree<Item_, Index_>::get_json_repr(json& json_repr, Itemizer itemizer) const
{
    json_repr["num_vertices"] = num_vertices();
    json_repr["num_levels"] = num_levels();

    std::vector<json> vertex_reprs;

    for (Index v = 0; v < num_vertices(); ++v)
    {
        vertex_reprs.emplace_back();
        json& vertex_repr = vertex_reprs.back();

        vertex_repr["vertex"] = v;
        vertex_repr["level"] = levels[v];
        vertex_repr["parent"] = parents[v];
        vertex_repr["children"] = children[v];
        itemizer(vertex_repr, vertices[v], *this);
    }

    json_repr["vertices"] = vertex_reprs;
}

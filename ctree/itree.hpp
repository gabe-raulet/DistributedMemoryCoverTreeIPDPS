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
    vertices.clear();
    levels.clear();
    children.clear();
    nlevels = 0;
}

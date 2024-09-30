template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::build(Real ghost_radius, Real split_ratio, Real switch_percent, Index min_hub_size, bool level_synch, bool threaded, bool verbose)
{
    using Hub = Hub<CoverTree>;
    using HubVector = typename Hub::HubVector;

    double t, elapsed;
    Real avg_hub_size;
    Index leaf_count, iter, num_hubs;
    HubVector hubs;

    Real switch_size = (switch_percent/100.0) * size;

    tree.clear();

    t = -omp_get_wtime();
    hubs.emplace_back(points);
    t += omp_get_wtime();
    elapsed = t;

    if (verbose)
    {
        fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] initialized root hub [farthest_point={},root_hub_radius={:.3f}]\n", __func__, elapsed, t, hubs.back().cand(), hubs.back().radius());
        std::cout << std::flush;
    }

    iter = 1;
    leaf_count = 0;
    std::vector<bool> leaf_flags(size, false);

    do
    {
        t = -omp_get_wtime();

        HubVector split_hubs, next_hubs;

        num_hubs = hubs.size();

        #pragma omp parallel if (threaded)
        {
            HubVector my_split_hubs, my_next_hubs;

            #pragma omp for nowait schedule(dynamic)
            for (Index i = 0; i < num_hubs; ++i)
            {
                Hub& hub = hubs[i];

                do { hub.add_new_leader(points); } while (level_synch && !hub.is_split(split_ratio));

                if (hub.is_split(split_ratio))
                {
                    hub.split_leaders(points);
                    hub.find_leaves(min_hub_size);
                    my_split_hubs.push_back(hub);
                }
                else
                {
                    my_next_hubs.push_back(hub);
                }
            }

            #pragma omp critical
            for (Hub& hub : my_next_hubs)
                next_hubs.push_back(hub);

            #pragma omp critical
            for (Hub& hub : my_split_hubs)
                split_hubs.push_back(hub);
        }

        for (Hub& split_hub : split_hubs)
        {
            leaf_count += split_hub.update_tree(tree, next_hubs, leaf_flags);
        }

        std::swap(hubs, next_hubs);

        t += omp_get_wtime();
        elapsed += t;

        avg_hub_size = (size-leaf_count+0.0)/hubs.size();

        if (verbose)
        {
            double leaf_percent = (100.0*leaf_count)/size;
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] {:.2f} percent leaves reached [iter={},levels={},vertices={},hubs={},avg_hub_size={:.3f}]\n", __func__, elapsed, t, leaf_percent, iter, tree.num_levels(), tree.num_vertices(), hubs.size(), avg_hub_size);
            std::cout << std::flush;
        }

        iter++;

    } while (static_cast<Real>(leaf_count) < switch_size && leaf_count < size);

    if (leaf_count < size) // need ghost trees
    {
        assert((!has_globids()));

        t = -omp_get_wtime();

        for (Hub& hub : hubs)
        {
            ghost_map[hub.repr()] = {ghost_map.size(), hub.add_hub_vertex(tree)};
            ghost_trees.emplace_back();
        }

        t += omp_get_wtime();
        elapsed += t;

        if (verbose)
        {
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] added {} remaining hub vertices to replication tree\n", __func__, elapsed, t, hubs.size());
            std::cout << std::flush;
        }

        t = -omp_get_wtime();


        #pragma omp parallel if (threaded)
        {
            IndexVectorVector my_ghost_hub_points(hubs.size());

            #pragma omp for nowait schedule(dynamic)
            for (Index i = 0; i < size; ++i)
            {
                if (leaf_flags[i]) continue;

                IndexVector hub_ids;
                hub_query(points[i], ghost_radius, hub_ids);

                for (Index hub_id : hub_ids)
                {
                    my_ghost_hub_points[ghost_map.at(hub_id).first].push_back(i);
                }
            }

            #pragma omp critical
            {
                for (Index i = 0; i < hubs.size(); ++i)
                {
                    for (Index id : my_ghost_hub_points[i])
                    {
                        ghost_trees[i].add_point(points[id], id);
                    }
                }
            }
        }

        t += omp_get_wtime();
        elapsed += t;

        if (verbose)
        {
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] setup {} ghost trees\n", __func__, elapsed, t, ghost_trees.size());
            std::cout << std::flush;
        }

        t = -omp_get_wtime();

        #pragma omp parallel for if (threaded) schedule(dynamic)
        for (Index i = 0; i < ghost_trees.size(); ++i)
        {
            ghost_trees[i].set_new_root(hubs[i].repr());
            ghost_trees[i].build(0.0, split_ratio, 100.0, min_hub_size, true, false, false);
        }

        t += omp_get_wtime();
        elapsed += t;

        if (verbose)
        {
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] built {} ghost trees\n", __func__, elapsed, t, ghost_trees.size());
            std::cout << std::flush;
        }
    }
}

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::point_query(const Point& query, Real epsilon, IndexVector& neighbors) const
{
    if (!has_ghost_trees()) reptree_point_query(query, epsilon, neighbors);
    else ghost_point_query(query, epsilon, neighbors);
}

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::ghost_point_query(const Point& query, Real epsilon, IndexVector& neighbors) const
{
    IndexVector hub_ids, es;

    hub_query(query, epsilon, hub_ids);
    reptree_point_query(query, epsilon, es);

    for (Index hub_id : hub_ids)
    {
        Index tree_slot = ghost_map.at(hub_id).first;
        ghost_trees[tree_slot].reptree_point_query(query, epsilon, es);
    }

    IndexSet allids(es.begin(), es.end());
    neighbors.assign(allids.begin(), allids.end());
}

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::reptree_point_query(const Point& query, Real epsilon, IndexVector& neighbors) const
{
    Index prev_size = neighbors.size();

    /*
     * Finds all points in tree that are in the `epsilon`-ball
     * centered at the query point.
     */

    const auto& [root_id, root_radius] = tree[0];

    /*
     * All descendants of the tree root are within the ball with
     * radius `root_radius` and center `root_pt`. If the query
     * point is not within a distance of `epsilon` of the boundary
     * of this ball, then it is impossible for their to be an
     * `epsilon` neighbor of the query in this tree.
     */
    if (distance(query, points[root_id]) > root_radius + epsilon)
        return;

    /*
     * Use a stack of tree vertices to implement recursive tree search.
     */
    IndexVector stack = {0};

    while (!stack.empty())
    {
        Index u = stack.back(); stack.pop_back();
        const auto& [uid, uradius] = tree[u];

        IndexVector children;
        tree.get_children(u, children);

        /*
         * If this vertex has no children then it is a leaf and
         * we check if it is an `epsilon`-neighbor of the query.
         */
        if (children.empty() && distance(query, points[uid]) <= epsilon)
        {
            neighbors.push_back(uid);
        }
        else
        {
            /*
             * Otherwise, we go through all the child vertices and
             * recursively search the sub-trees that could have `epsilon`-neighbors
             * as descendants.
             */
            for (Index v : children)
            {
                const auto& [vid, vradius] = tree[v];

                if (distance(query, points[vid]) <= vradius + epsilon)
                    stack.push_back(v);
            }
        }
    }

    if (has_globids()) std::for_each(neighbors.begin() + prev_size, neighbors.end(), [&](Index& id) { id = globids[id]; });
}

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::hub_query(const Point& query, Real ghost_radius, IndexVector& hub_ids) const
{
    IndexVector stack = {0};

    while (!stack.empty())
    {
        Index u = stack.back(); stack.pop_back();
        const auto& [uid, uradius] = tree[u];

        IndexVector children;
        tree.get_children(u, children);

        auto it = ghost_map.find(uid);

        if (it != ghost_map.end())
        {
            const auto& [_, vtx] = it->second;

            if (u == vtx && distance(query, points[uid]) <= uradius + ghost_radius)
                hub_ids.push_back(uid);
        }


        for (Index v : children)
        {
            const auto& [vid, vradius] = tree[v];

            if (distance(query, points[vid]) <= vradius + ghost_radius)
                stack.push_back(v);
        }
    }
}

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::write_tree_file(const char *fname) const
{
    FILE *f = fopen(fname, "w");

    Index nlevels = tree.num_levels();
    Index nverts = tree.num_vertices();
    Index npoints = num_points();

    fmt::print(f, "vertex_id\tparent_id\tpoint_id\tradius\n");
    fflush(f);

    for (Index i = 0; i < nverts; ++i)
    {
        Index parent_id = tree.parents[i];
        Index point_id = tree.vertices[i].id;
        Real radius = tree.vertices[i].radius;

        fmt::print(f, "{}\t{}\t{}\t{:.5f}\n", i, parent_id, point_id, radius);
        fflush(f);
    }

    fclose(f);
}

template <class PointTraits_, class Distance_, index_type Index_>
bool CoverTree<PointTraits_, Distance_, Index_>::is_correct(Real split_ratio) const
{
    const Index n = tree.num_vertices();
    std::vector<bool> ptflags(size, false);

    Index found = 0;

    /*
     * Go through each vertex in tree
     */
    for (Index u = 0; u < n; ++u)
    {
        IndexVector children;
        tree.get_children(u, children);

        const auto& [uid, urad] = tree[u];

        if (!ptflags[uid]) ptflags[uid] = true, found++;

        for (Index v : children)
        {
            const auto& [vid, vrad] = tree[v];

            /*
             * Child point must be within the ball of the parent
             * or we fail the covering condition
             */

            if (distance(points[uid], points[vid]) > urad)
            {
                fmt::print(stderr, "[err::{}] failed covering condition!\n", __func__);
                return false;
            }
        }

        /*
         * Pairwise compare distinct children vertices (siblings)
         */
        for (Index p1 : children)
            for (Index p2 : children)
            {
                if (p1 == p2 || tree.is_leaf(p1) || tree.is_leaf(p2))
                    continue;

                const auto& [id1, rad1] = tree[p1];
                const auto& [id2, rad2] = tree[p2];

                /*
                 * Sibling points must be separated by a constant
                 * ratio (split_ratio) of their parent's ball radius
                 */
                if (distance(points[id1], points[id2]) <= urad*split_ratio)
                {
                    fmt::print(stderr, "[err::{}] failed sibling separation condition!\n", __func__);
                    return false;
                }
            }
    }

    if (found != size)
    {
        fmt::print(stderr, "[err::{}] failed covering condition\n", __func__);
        return false;
    }

    return true;
}

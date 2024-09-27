template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::build(const PointVector& points, Real split_ratio, Real switch_size, Index min_hub_size, bool level_synch, bool threaded, bool verbose)
{
    using Hub = Hub<CoverTree>;
    using HubVector = typename Hub::HubVector;

    double t, elapsed;
    Real avg_hub_size;
    BallTree balltree;
    Index leaf_count, iter, num_hubs;
    IndexVector hub_map;
    HubVector hubs;

    tree.clear();
    size = points.size();

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
    hub_map.resize(size, -1);

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
            leaf_count += split_hub.update_tree(balltree, next_hubs, points, hub_map);
        }

        std::swap(hubs, next_hubs);

        t += omp_get_wtime();
        elapsed += t;

        avg_hub_size = (size-leaf_count+0.0)/hubs.size();

        if (verbose)
        {
            double leaf_percent = (100.0*leaf_count)/size;
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] {:.2f} percent leaves reached on iteration {} [vertices={},levels={},hubs={},avg_hub_size={:.3f}]\n", __func__, elapsed, t, leaf_percent, iter, balltree.num_vertices(), balltree.num_levels(), hubs.size(), avg_hub_size);
            std::cout << std::flush;
        }

        iter++;

    } while (avg_hub_size > switch_size && leaf_count < size);

    t = -omp_get_wtime();
    fill_point_ball_tree(balltree, points);
    t += omp_get_wtime();
    elapsed += t;

    if (verbose)
    {
        fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] added points to tree\n", __func__, elapsed, t);
        std::cout << std::flush;
    }

}

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::point_query(const Point& query, Real epsilon, IndexVector& neighbors) const
{
    /*
     * Finds all points in tree that are in the `epsilon`-ball
     * centered at the query point.
     */

    const auto& [root_pt, root_id, root_radius] = tree[0];

    /*
     * All descendants of the tree root are within the ball with
     * radius `root_radius` and center `root_pt`. If the query
     * point is not within a distance of `epsilon` of the boundary
     * of this ball, then it is impossible for their to be an
     * `epsilon` neighbor of the query in this tree.
     */
    if (distance(query, root_pt) > root_radius + epsilon)
        return;

    /*
     * Use a stack of tree vertices to implement recursive tree search.
     */
    IndexVector stack = {0};

    while (!stack.empty())
    {
        Index u = stack.back(); stack.pop_back();
        const auto& [upt, uid, uradius] = tree[u];

        IndexVector children;
        tree.get_children(u, children);

        /*
         * If this vertex has no children then it is a leaf and
         * we check if it is an `epsilon`-neighbor of the query.
         */
        if (children.empty() && distance(query, upt) <= epsilon)
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
                const auto& [vpt, vid, vradius] = tree[v];

                if (distance(query, vpt) <= vradius + epsilon)
                    stack.push_back(v);
            }
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

        const auto& [upt, uid, urad] = tree[u];

        if (!ptflags[uid]) ptflags[uid] = true, found++;

        for (Index v : children)
        {
            const auto& [vpt, vid, vrad] = tree[v];

            /*
             * Child point must be within the ball of the parent
             * or we fail the covering condition
             */

            if (distance(upt, vpt) > urad)
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

                const auto& [pt1, id1, rad1] = tree[p1];
                const auto& [pt2, id2, rad2] = tree[p2];

                /*
                 * Sibling points must be separated by a constant
                 * ratio (split_ratio) of their parent's ball radius
                 */
                if (distance(pt1, pt2) <= urad*split_ratio)
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

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::fill_point_ball_tree(const BallTree& balltree, const PointVector& points)
{
    tree.levels = std::move(balltree.levels);
    tree.parents = std::move(balltree.parents);
    tree.children = std::move(balltree.children);
    tree.nlevels = balltree.nlevels;

    Index n = balltree.num_vertices();
    tree.vertices.resize(n);

    #pragma omp parallel for
    for (Index i = 0; i < n; ++i)
    {
        const Ball& ball = balltree[i];
        tree.vertices[i] = {points[ball.id], ball.id, ball.radius};
    }
}

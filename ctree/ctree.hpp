template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::build(Real ghost_radius, Real split_ratio, Real switch_percent, Index min_hub_size, bool level_synch, bool threaded, bool verbose)
{
    using Hub = Hub<CoverTree>;
    using HubVector = typename Hub::HubVector;

    double t, elapsed;
    Real avg_hub_size;
    BallTree balltree;
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

    IndexVector pt2hub(size, 0);
    IndexVector hub2vtx(size, -1);
    std::vector<bool> is_leaf(size, false);

#ifdef LOG
    auto itemizer = [&] (json& vertex_repr, const Ball& ball, const BallTree& balltree)
    {
        vertex_repr["point"] = ball.id;
        vertex_repr["radius"] = ball.radius;
        vertex_repr["is_leaf"] = vertex_repr["children"].empty() && is_leaf[ball.id];
        vertex_repr["hub"] = pt2hub[ball.id];
        vertex_repr["nested"] = balltree.vertices[vertex_repr["parent"]].id == ball.id;
    };

    std::vector<json> iterations;

#endif

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
            hub2vtx[split_hub.repr()] = split_hub.add_hub_vertex(balltree);

            for (const auto& hub_point : split_hub.get_hub_points())
                pt2hub[hub_point.id] = split_hub.repr();
        }

        for (Hub& split_hub : split_hubs)
        {
            Index num_leaves = split_hub.add_hub_leaves(balltree, next_hubs);
            leaf_count += num_leaves;

            for (Index leaf : split_hub.get_leaves())
                is_leaf[leaf] = true;
        }

        std::swap(hubs, next_hubs);

        t += omp_get_wtime();
        elapsed += t;

        avg_hub_size = (size-leaf_count+0.0)/hubs.size();

        if (verbose)
        {
            double leaf_percent = (100.0*leaf_count)/size;
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] {:.2f} percent leaves reached [iter={},levels={},vertices={},hubs={},avg_hub_size={:.3f}]\n", __func__, elapsed, t, leaf_percent, iter, balltree.num_levels(), balltree.num_vertices(), hubs.size(), avg_hub_size);
            std::cout << std::flush;
        }

#ifdef LOG

        iterations.emplace_back();
        json& iter_json = iterations.back();
        iter_json["iter"] = iter;
        balltree.get_json_repr(iter_json["tree"], itemizer);
#endif

        iter++;

    } while (static_cast<Real>(leaf_count) < switch_size && leaf_count < size);

    t = -omp_get_wtime();
    fill_point_ball_tree(balltree, points);
    t += omp_get_wtime();
    elapsed += t;

    if (verbose)
    {
        fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] added points to tree\n", __func__, elapsed, t);
        std::cout << std::flush;
    }

#ifdef LOG

    if (!has_globids())
    {
        json tree_repr;
        tree_repr["iterations"] = iterations;
        std::ofstream os("tree_repr.json");
        os << std::setw(4) << tree_repr << std::endl;
        os.close();
    }

#endif

    if (leaf_count < size) // need ghost trees
    {
        assert((!has_globids()));

        IndexMap hub_slots;

        t = -omp_get_wtime();

        for (const Hub& hub : hubs)
        {
            hub_slots[hub.repr()] = hub_slots.size();

            Index relroot = -1;
            PointVector hub_points;
            IndexVector hub_point_ids;
            hub_point_ids.reserve(hub.size());

            for (const auto& hub_point : hub.get_hub_points())
            {
                if (hub_point.id == hub.repr())
                {
                    relroot = hub_point_ids.size();
                }

                hub_point_ids.push_back(hub_point.id);
                hub_points.push_back(points[hub_point.id]);
            }

            assert((relroot >= 0 && relroot < hub.size()));

            if (relroot != 0)
            {
                std::swap(hub_point_ids.front(), hub_point_ids[relroot]);
                std::swap(hub_points.front(), hub_points[relroot]);
            }

            ghost_trees.emplace_back(hub_points, hub_point_ids);
        }

        t += omp_get_wtime();
        elapsed += t;

        if (verbose)
        {
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] setup {} ghost trees\n", __func__, elapsed, t, ghost_trees.size());
            std::cout << std::flush;
        }

        t = -omp_get_wtime();

        #pragma omp parallel for
        for (Index i = 0; i < ghost_trees.size(); ++i)
        {
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

    if (has_globids())
    {
        std::for_each(neighbors.begin(), neighbors.end(), [&](Index& id) { id = globids[id]; });
    }
    else if (has_ghost_trees()) // this and has_globids() should be mutually exclusive
    {
        //IndexVector ghost_neighbors;
        //IndexSet all_neighbors(neighbors.begin(), neighbors.end());

        //for (const auto& ghost_tree : ghost_trees)
        //{
        //    ghost_neighbors.clear();
        //    ghost_tree.point_query(query, epsilon, ghost_neighbors);
        //    all_neighbors.insert(ghost_neighbors.begin(), ghost_neighbors.end());
        //}

        //neighbors.assign(all_neighbors.begin(), all_neighbors.end());
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

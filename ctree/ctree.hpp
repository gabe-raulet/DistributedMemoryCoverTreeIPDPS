template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::build(Real split_ratio, Index min_hub_size, bool threaded, bool verbose)
{
    using Hub = Hub<CoverTree>;
    using HubVector = typename Hub::HubVector;

    Index size = num_points();

    double t, elapsed;
    Real avg_hub_size;
    Index leaf_count, iter, num_hubs;
    HubVector hubs;

    BallTree balltree;

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

        HubVector next_hubs;

        num_hubs = hubs.size();

        #pragma omp parallel if (threaded)
        {
            HubVector my_hubs;

            #pragma omp for nowait schedule(dynamic)
            for (Index i = 0; i < num_hubs; ++i)
            {
                Hub& hub = hubs[i];

                do { hub.add_new_leader(points); } while (!hub.is_split(split_ratio));

                hub.split_leaders(points);
                hub.find_leaves(min_hub_size);
                my_hubs.push_back(hub);
            }

            #pragma omp critical
            for (Hub& hub : my_hubs)
            {
                leaf_count += hub.update_tree(balltree, next_hubs, leaf_flags);
            }
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

        iter++;

    } while (leaf_count < size);

    t = -omp_get_wtime();

    auto itemizer = [&](const Ball& ball) -> PointBall { return {points[ball.id], ball.id, ball.radius}; };

    balltree.itemize_new_tree(tree, itemizer);

    t += omp_get_wtime();
    elapsed += t;

    if (verbose)
    {
        fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] added points to tree\n", __func__, elapsed, t);
        std::cout << std::flush;
    }
}

template <class PointTraits_, class Distance_, index_type Index_>
bool CoverTree<PointTraits_, Distance_, Index_>::is_correct(Real split_ratio) const
{
    const Index n = tree.num_vertices();
    std::vector<bool> ptflags(num_points(), false);

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

    if (found != num_points())
    {
        fmt::print(stderr, "[err::{}] failed covering condition\n", __func__);
        return false;
    }

    return true;
}

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::point_query(const Point& query, Real epsilon, IndexVector& neighbors) const
{
    Index prev_size = neighbors.size();

    const auto& [root_pt, root_id, root_radius] = tree[0];

    if (distance(query, root_pt) > root_radius + epsilon)
        return;

    IndexVector stack = {0};

    while (!stack.empty())
    {
        Index u = stack.back(); stack.pop_back();
        const auto& [upt, uid, uradius] = tree[u];

        IndexVector children;
        tree.get_children(u, children);

        if (children.empty() && distance(query, upt) <= epsilon)
        {
            neighbors.push_back(uid);
        }
        else
        {
            for (Index v : children)
            {
                const auto& [vpt, vid, vradius] = tree[v];

                if (distance(query, vpt) <= vradius + epsilon)
                    stack.push_back(v);
            }
        }
    }

    if (has_globids()) std::for_each(neighbors.begin()+prev_size, neighbors.end(), [&](Index& id) { id = globids[id]; });
}

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::add_point(Point pt, Index globid)
{
    points.push_back(pt);
    globids.push_back(globid);
}

template <class PointTraits_, class Distance_, index_type Index_>
void CoverTree<PointTraits_, Distance_, Index_>::set_new_root(Index root)
{
    if (has_globids())
    {
        Index where = std::find(globids.begin(), globids.end(), root) - globids.begin();

        if (where != 0)
        {
            std::swap(globids[0], globids[where]);
            std::swap(points[0], points[where]);
        }
    }
    else if (root != 0) std::swap(points[0], points[root]);
}

template <class PointTraits_, class Distance_, index_type Index_>
typename CoverTree<PointTraits_, Distance_, Index_>::Index
CoverTree<PointTraits_, Distance_, Index_>::build_epsilon_graph(Real radius, IndexVectorVector& neighbors) const
{
    Index size = num_points();
    neighbors.resize(size, {});

    Index num_edges = 0;

    #pragma omp parallel for reduction(+:num_edges)
    for (Index i = 0; i < size; ++i)
    {
        point_query(points[i], radius, neighbors[i]);
        num_edges += neighbors[i].size();
    }

    return num_edges;
}


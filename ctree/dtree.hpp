template <class PointTraits_, class Distance_, index_type Index_>
DistCoverTree<PointTraits_, Distance_, Index_>::DistCoverTree(const PointVector& mypoints, const Comm& comm) : comm(comm), mysize(mypoints.size()), mypoints(mypoints), offsets(comm.size())
{
    IndexVector counts(comm.size());
    counts[comm.rank()] = mysize;
    comm.allgather(counts);

    std::exclusive_scan(counts.begin(), counts.end(), offsets.begin(), static_cast<Index>(0));

    myoffset = offsets[comm.rank()];
    totsize = offsets.back() + counts.back();
}

template <class PointTraits_, class Distance_, index_type Index_>
void DistCoverTree<PointTraits_, Distance_, Index_>::build(Real radius, Real split_ratio, Real switch_percent, Index min_hub_size, bool verbose)
{
    using DistHub = DistHub<DistCoverTree>;
    using DistHubVector = typename DistHub::DistHubVector;

    auto timer = comm.get_timer();
    double elapsed,  t;
    Real avg_hub_size;
    Index leaf_count, iter, num_hubs;
    DistHubVector hubs;
    BallTree repballtree;

    Real switch_size = (switch_percent/100.0) * totsize;

    timer.start_timer();

    DistHub::init_mpi_argmax_op();

    Point root_pt = mypoints.front();
    comm.bcast(root_pt, 0);

    hubs.emplace_back(mypoints, root_pt, *this);

    timer.stop_timer();
    t = elapsed = timer.get_max_time();

    if (verbose && !comm.rank())
    {
        fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] initialized root hub [farthest_point={},root_hub_radius={:.3f}]\n", __func__, elapsed, t, hubs.back().cand(), hubs.back().radius());
        std::cout << std::flush;
    }

    iter = 1;
    leaf_count = 0;

    do
    {
        timer.start_timer();

        DistHubVector split_hubs, next_hubs;
        num_hubs = hubs.size();

        for (Index i = 0; i < num_hubs; ++i)
            hubs[i].add_new_leader(mypoints);

        DistHub::synchronize_hubs(hubs, comm);

        for (DistHub& hub : hubs)
        {
            if (hub.is_split(split_ratio))
                split_hubs.push_back(hub);
            else
                next_hubs.push_back(hub);
        }

        for (DistHub& hub : split_hubs)
        {
            hub.split_leaders(mypoints);
        }

        DistHub::synchronize_split_hubs(split_hubs, min_hub_size, comm);

        for (DistHub& split_hub : split_hubs)
        {
            leaf_count += split_hub.update_tree(repballtree, next_hubs);
        }

        std::swap(hubs, next_hubs);

        timer.stop_timer();
        t = timer.get_max_time();
        elapsed += t;

        avg_hub_size = (totsize - leaf_count + 0.0) / hubs.size();

        if (verbose && !comm.rank())
        {
            double leaf_percent = (100.0*leaf_count)/totsize;
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] {:.2f} percent leaves reached [iter={},levels={},vertices={},hubs={},avg_hub_size={:.3f}]\n", __func__, elapsed, t, leaf_percent, iter, repballtree.num_levels(), repballtree.num_vertices(), hubs.size(), avg_hub_size);
            std::cout << std::flush;
        }

        iter++;

    } while (static_cast<Real>(leaf_count) < switch_size && leaf_count < totsize);

    assert((leaf_count <= totsize));

    using PointTriple = std::tuple<Index, Index, Point>; // hub id, point id, point
    using PointTripleVector = std::vector<PointTriple>;

    std::vector<PointTripleVector> sendbufs(comm.size());
    PointTripleVector recvbuf;

    if (leaf_count < totsize)
    {
        timer.start_timer();

        for (DistHub& hub : hubs)
        {
            Index repr = hub.repr();
            int owner = point_owner(repr);

            if (owner == comm.rank())
            {
                ghost_trees.try_emplace(repr);
            }

            ghost_map[repr] = hub.add_hub_vertex(repballtree);
            hub_sizes[repr] = hub.size();

            for (const auto& hub_point : hub.get_hub_points())
            {
                Index p = hub_point.id;
                Point pt = mypoints[p-myoffset];
                sendbufs[owner].emplace_back(repr, p, pt);
            }
        }

        timer.stop_timer();
        t = timer.get_max_time();
        elapsed += t;

        if (verbose && !comm.rank())
        {
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] initialized {} local trees and loaded outgoing hub point buffers\n", __func__, elapsed, t, hubs.size());
            std::cout << std::flush;
        }
    }

    timer.start_timer();

    PointPairVector my_point_pairs, point_pairs;

    for (Index i = 0; i < repballtree.num_vertices(); ++i)
    {
        Index p = repballtree[i].id;

        if (repballtree.is_leaf(i) && owns_point(p))
        {
            my_point_pairs.emplace_back(p, mypoints[p-myoffset]);
        }
    }

    comm.allgatherv(my_point_pairs, point_pairs);

    PointMap reptree_points;
    reptree_points.reserve(point_pairs.size());

    for (const auto& [globid, pt] : point_pairs)
        reptree_points.insert({globid, pt});

    auto itemizer = [&](const Ball& ball) -> PointBall { return {reptree_points.at(ball.id), ball.id, ball.radius}; };

    repballtree.itemize_new_tree(reptree, itemizer);

    timer.stop_timer();
    t = timer.get_max_time();
    elapsed += t;

    if (verbose && !comm.rank())
    {
        fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] added points to replication tree vertices\n", __func__, elapsed, t);
        std::cout << std::flush;
    }

    if (leaf_count == totsize)
    {
        DistHub::free_mpi_argmax_op();
        return;
    }

    timer.start_timer();

    comm.alltoallv(sendbufs, recvbuf);

    /* std::sort(recvbuf.begin(), recvbuf.end()); */

    for (const auto& [repr, p, pt] : recvbuf)
    {
        /* assert((ghost_trees.contains(repr))); */
        ghost_trees[repr].add_point(pt, p);
    }

    timer.stop_timer();
    t = timer.get_max_time();
    elapsed += t;


    if (verbose && !comm.rank())
    {
        fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] collected hubs to local ranks via alltoall\n", __func__, elapsed, t);
        std::cout << std::flush;
    }

    double ranktime, roottime;

    timer.start_timer();
    ranktime = -MPI_Wtime();

    recvbuf.clear();
    std::for_each(sendbufs.begin(), sendbufs.end(), [](auto& sendbuf) { sendbuf.clear(); });

    Index my_num_queries = 0, root_num_queries = 0;

    for (auto& [repr, ghost_tree] : ghost_trees)
    {
        /* ghost_tree.set_new_root(repr); */
        Index n = ghost_tree.num_points();
        const Point* tree_points = ghost_tree.point_data();
        const Index* tree_globids = ghost_tree.globid_data();

        for (Index i = 0; i < n; ++i)
        {
            Point pt = tree_points[i];
            Index p = tree_globids[i];
            my_num_queries++;

            IndexVector neighbors, ghost_hubs;
            point_query(pt, radius, neighbors, ghost_hubs, repr);

            for (Index ghost_hub : ghost_hubs)
            {
                sendbufs[point_owner(ghost_hub)].emplace_back(ghost_hub, p, pt);
            }
        }
    }

    ranktime += MPI_Wtime();

    if (comm.rank()==comm.size()-1)
    {
        roottime = -MPI_Wtime();

        for (const auto& [p, pt] : reptree_points)
        {
            if (ghost_map.contains(p))
                continue;

            root_num_queries++;

            IndexVector neighbors, ghost_hubs;
            point_query(pt, radius, neighbors, ghost_hubs);

            for (Index ghost_hub : ghost_hubs)
            {
                sendbufs[point_owner(ghost_hub)].emplace_back(ghost_hub, p, pt);
            }
        }

        roottime += MPI_Wtime();
    }

    if (verbose)
    {
        fmt::print("[msg::{},elapsed={:.3f},ranktime={:.3f}] rank {} finished querying {} hub points against replication tree\n", __func__, ranktime+elapsed, ranktime, comm.rank(), my_num_queries);
        if (comm.rank()==comm.size()-1) fmt::print("[msg::{},elapsed={:.3f},ranktime={:.3f}] rank {} finished querying {} replication-tree leaves against replication tree\n", __func__, ranktime+roottime+elapsed, roottime, comm.rank(), root_num_queries);
        std::cout << std::flush;
    }

    comm.alltoallv(sendbufs, recvbuf);

    /* std::sort(recvbuf.begin(), recvbuf.end()); */

    for (const auto& [repr, p, pt] : recvbuf)
    {
        /* assert((ghost_trees.contains(repr))); */
        ghost_trees[repr].add_point(pt, p);
    }

    timer.stop_timer();
    t = timer.get_max_time();
    elapsed += t;

    if (verbose && !comm.rank())
    {
        fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] finished setting up {} ghost trees\n", __func__, elapsed, t, hubs.size());
        std::cout << std::flush;
    }

    timer.start_timer();
    ranktime = -MPI_Wtime();

    for (auto& [repr, ghost_tree] : ghost_trees)
    {
        ghost_tree.build(split_ratio, min_hub_size, false, false);
    }

    ranktime += MPI_Wtime();

    if (verbose)
    {
        fmt::print("[msg::{},elapsed={:.3f},ranktime={:.3f}] rank {} finished building {} ghost trees\n", __func__, ranktime+elapsed, ranktime, comm.rank(), ghost_trees.size());
        std::cout << std::flush;
    }

    DistHub::free_mpi_argmax_op();
}


template <class PointTraits_, class Distance_, index_type Index_>
typename DistCoverTree<PointTraits_, Distance_, Index_>::Index
DistCoverTree<PointTraits_, Distance_, Index_>::build_epsilon_graph(Real radius, IndexVectorVector& myneighbors) const
{
    myneighbors.resize(mysize, {});
    Index num_edges = 0;
    comm.allreduce(num_edges, MPI_SUM);
    return num_edges;
}

template <class PointTraits_, class Distance_, index_type Index_>
void DistCoverTree<PointTraits_, Distance_, Index_>::point_query(const Point& query, Real radius, IndexVector& neighbors, IndexVector& ghost_hubs, Index query_hub) const
{
    IndexVector stack = {0};

    while (!stack.empty())
    {
        Index u = stack.back(); stack.pop_back();
        const auto& [upt, uid, uradius] = reptree[u];

        IndexVector children;
        reptree.get_children(u, children);

        if (!children.empty())
        {
            for (Index v : children)
            {
                const auto& [vpt, vid, vradius] = reptree[v];

                if (distance(query, vpt) <= vradius + radius)
                    stack.push_back(v);
            }
        }
        else
        {
            Real dist = distance(query, upt);

            if (uid != query_hub && ghost_map.contains(uid) && u == ghost_map.at(uid) && dist <= uradius + radius)
            {
                ghost_hubs.push_back(uid);
            }
            else if (dist <= radius)
            {
                neighbors.push_back(uid);
            }
        }
    }
}

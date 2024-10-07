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

        if (repballtree.is_leaf(i))
        {
            if (!ghost_map.contains(p))
                rep_leaves.push_back(p);

            if (owns_point(p))
            {
                my_point_pairs.emplace_back(p, mypoints[p-myoffset]);
            }
        }
    }

    comm.allgatherv(my_point_pairs, point_pairs);

    reptree_points.clear();
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

    Index my_num_queries = 0, rep_num_queries = 0;

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

    IndexVector rep_counts(comm.size()), rep_offsets(comm.size());
    get_balanced_counts(rep_counts, point_pairs.size());

    std::exclusive_scan(rep_counts.begin(), rep_counts.end(), rep_offsets.begin(), static_cast<Index>(0));

    for (Index i = 0; i < rep_counts[comm.rank()]; ++i)
    {
        const auto& [p, pt] = point_pairs[i+rep_offsets[comm.rank()]];

        if (ghost_map.contains(p))
            continue;

        rep_num_queries++;

        IndexVector neighbors, ghost_hubs;
        point_query(pt, radius, neighbors, ghost_hubs);

        for (Index ghost_hub : ghost_hubs)
        {
            sendbufs[point_owner(ghost_hub)].emplace_back(ghost_hub, p, pt);
        }
    }

    ranktime += MPI_Wtime();

    if (verbose)
    {
        fmt::print("[msg::{},elapsed={:.3f},ranktime={:.3f}] rank {} finished querying {} hub points and {} replication-leaf points against replication tree\n", __func__, ranktime+elapsed, ranktime, comm.rank(), my_num_queries, rep_num_queries);
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

    Index num_hub_points = 0, num_tree_points = 0;

    for (auto& [repr, ghost_tree] : ghost_trees)
    {
        ghost_tree.build(split_ratio, min_hub_size, false, false);
        num_hub_points += hub_sizes.at(repr);
        num_tree_points += ghost_tree.num_points();
    }

    ranktime += MPI_Wtime();

    if (verbose)
    {
        fmt::print("[msg::{},elapsed={:.3f},ranktime={:.3f}] rank {} finished building {} ghost trees [num_hub_points={},num_tree_points={},percent_ghost_points={:.3f}]\n", __func__, ranktime+elapsed, ranktime, comm.rank(), ghost_trees.size(), num_hub_points, num_tree_points, 100.*(num_tree_points-num_hub_points) / num_tree_points);
        std::cout << std::flush;
    }

    DistHub::free_mpi_argmax_op();
}


template <class PointTraits_, class Distance_, index_type Index_>
typename DistCoverTree<PointTraits_, Distance_, Index_>::Index
DistCoverTree<PointTraits_, Distance_, Index_>::build_epsilon_graph(Real radius, IndexVectorVector& myneighbors) const
{
    auto timer = comm.get_timer();
    timer.start_timer();

    int myrank = comm.rank();
    int nprocs = comm.size();

    myneighbors.resize(mysize, {});
    Index num_edges = 0;

    using Edge = std::pair<Index, Index>;
    using EdgeVector = std::vector<Edge>;

    std::vector<EdgeVector> sendbufs(nprocs);
    EdgeVector recvbuf;

    IndexSet rep_leaves_set(rep_leaves.begin(), rep_leaves.end());
    IndexVector rep_leaves_counts(nprocs), rep_leaves_offsets(nprocs);

    get_balanced_counts(rep_leaves_counts, rep_leaves.size());
    std::exclusive_scan(rep_leaves_counts.begin(), rep_leaves_counts.end(), rep_leaves_offsets.begin(), static_cast<Index>(0));

    for (const auto& [repr, ghost_tree] : ghost_trees)
    {
        Index n = hub_sizes.at(repr);
        const Point* pts = ghost_tree.point_data();
        const Index* ids = ghost_tree.globid_data();

        for (Index i = 0; i < n; ++i)
        {
            Point pt = pts[i];
            Index p = ids[i];

            IndexVector neighbors;
            ghost_tree.point_query(pt, radius, neighbors);

            int owner = point_owner(p);

            for (Index dest : neighbors)
            {
                sendbufs[owner].emplace_back(p, dest);

                if (rep_leaves_set.contains(dest))
                {
                    sendbufs[point_owner(dest)].emplace_back(dest, p);
                }
            }
        }
    }

    auto it = rep_leaves.begin() + rep_leaves_offsets[myrank];

    for (Index i = 0; i < rep_leaves_counts[myrank]; ++i)
    {
        Index p = *it++;
        Point pt = reptree_points.at(p);
        int owner = point_owner(p);

        IndexVector neighbors, dummy;
        point_query(pt, radius, neighbors, dummy);

        for (Index dest : neighbors)
        {
            sendbufs[owner].emplace_back(p, dest);
        }
    }

    timer.stop_timer();
    double t = timer.get_max_time();

    if (!comm.rank())
    {
        fmt::print("[msg::{},time={:.3f}] finished finding local neighbors\n", __func__, t);
        std::cout << std::flush;
    }

    timer.start_timer();
    comm.alltoallv(sendbufs, recvbuf);
    timer.stop_timer();

    t = timer.get_max_time();

    if (!comm.rank())
    {
        fmt::print("[msg::{},time={:.3f}] finished sending neighbors to their owners\n", __func__, t);
        std::cout << std::flush;
    }

    timer.start_timer();
    for (const auto& [u, v] : recvbuf)
    {
        myneighbors[u-myoffset].push_back(v);
        num_edges++;
    }

    comm.allreduce(num_edges, MPI_SUM);
    timer.stop_timer();
    t = timer.get_max_time();

    if (!comm.rank())
    {
        fmt::print("[msg::{},time={:.3f}] finished receiving neighbors\n", __func__, t);
        std::cout << std::flush;
    }

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

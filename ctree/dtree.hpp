template <class PointTraits_, class Distance_, index_type Index_>
DistCoverTree<PointTraits_, Distance_, Index_>::DistCoverTree(const PointVector& mypoints, const Comm& comm) : comm(comm), mysize(mypoints.size()), mypoints(mypoints)
{
    comm.exscan(mysize, myoffset, MPI_SUM, (Index)0);

    totsize = myoffset + mysize;
    comm.bcast(totsize, comm.size()-1);
}

template <class PointTraits_, class Distance_, index_type Index_>
DistCoverTree<PointTraits_, Distance_, Index_>::DistCoverTree(const PointVector& points, int root, const Comm& comm) : comm(comm)
{
    std::vector<int> sendcounts;

    if (comm.rank() == root)
    {
        sendcounts.resize(comm.size());
        get_balanced_counts(sendcounts, (size_t)points.size());
    }

    comm.scatterv(points, sendcounts, mypoints, root);
    mysize = mypoints.size();

    comm.exscan(mysize, myoffset, MPI_SUM, (Index)0);

    totsize = myoffset + mysize;
    comm.bcast(totsize, comm.size()-1);
}

template <class PointTraits_, class Distance_, index_type Index_>
void DistCoverTree<PointTraits_, Distance_, Index_>::build(Real ghost_radius, Real split_ratio, Real switch_percent, Index min_hub_size, bool verbose)
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
    IndexSet leaf_pts;

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
            leaf_count += split_hub.update_tree(repballtree, next_hubs, leaf_pts);
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

    if (leaf_count < totsize) // need ghost trees
    {
        timer.start_timer();

        for (DistHub& hub : hubs)
        {
            ghost_map[hub.repr()] = {ghost_map.size(), hub.add_hub_vertex(repballtree)};
        }

        timer.stop_timer();
        t = timer.get_max_time();
        elapsed += t;

        if (verbose && !comm.rank())
        {
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] added {} remaining hub vertices to skeleton tree\n", __func__, elapsed, t, hubs.size());
            std::cout << std::flush;
        }
    }

    timer.start_timer();

    Index num_rep_points = build_replication_tree(repballtree);

    timer.stop_timer();
    t = timer.get_max_time();
    elapsed += t;

    if (verbose && !comm.rank())
    {
        fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] built replication tree from skeleton tree [rep_vertices={},rep_levels={},rep_points={},rep_avg_nesting={:.3f}]\n", __func__, elapsed, t, reptree.num_vertices(), reptree.num_levels(), num_rep_points, (reptree.num_vertices()+0.0)/num_rep_points);
        std::cout << std::flush;
    }

    if (leaf_count < totsize) // need ghost trees
    {
        timer.start_timer();

        IndexVectorVector my_ghost_hub_points(hubs.size());

        for (Index i = 0; i < mysize; ++i)
        {
            if (leaf_pts.contains(myoffset+i))
                continue;

            IndexVector hub_ids;
            hub_query(mypoints[i], ghost_radius, hub_ids);

            for (Index hub_id : hub_ids)
            {
                my_ghost_hub_points[ghost_map.at(hub_id).first].push_back(myoffset+i);
            }
        }

        timer.stop_timer();
        t = timer.get_max_time();
        elapsed += t;

        if (verbose && !comm.rank())
        {
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] queried point partitions against replication tree\n", __func__, elapsed, t);
            std::cout << std::flush;
        }

        timer.start_timer();

        IndexVector offsets(comm.size());
        offsets[comm.rank()] = myoffset;
        comm.allgather(offsets);

        IndexMap hub_to_proc_map; /* (global) maps hub representatives to their processor owners */

         // go through each hub
        for (const auto& hub : hubs)
        {
            Index repr = hub.repr();
            Index where = (std::upper_bound(offsets.begin(), offsets.end(), repr) - offsets.begin()) - 1;
            hub_to_proc_map.insert({repr, where}); // map the hub to its destination processor

            if (myoffset <= repr && repr < myoffset + mysize)
                ghost_trees.try_emplace(repr); // initialize local hub ghost tree
        }

        using PointTriple = std::tuple<Index, Index, Point>; /* point id, hub repr, point */
        using PointTripleVector = std::vector<PointTriple>;

        std::vector<PointTripleVector> sendbufs(comm.size());
        PointTripleVector recvbuf;

        // go through each hub
        for (Index slot = 0; slot < hubs.size(); ++slot)
        {
            // go through each local hub point that belongs to the augmented ghost hub
            for (Index id : my_ghost_hub_points[slot])
            {
                // determine where that hub is supposed to belong
                int where = hub_to_proc_map.at(hubs[slot].repr());
                sendbufs[where].emplace_back(id, hubs[slot].repr(), mypoints[id-myoffset]); // send point to its destination
            }
        }

        comm.alltoallv(sendbufs, recvbuf);

        for (const auto& [id, repr, pt] : recvbuf)
        {
            ghost_trees[repr].add_point(pt, id);
        }

        timer.stop_timer();
        t = timer.get_max_time();
        elapsed += t;

        if (verbose && !comm.rank())
        {
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] constructed local ghost hubs via alltoall\n", __func__, elapsed, t);
            std::cout << std::flush;
        }
    }

    DistHub::free_mpi_argmax_op();
}

template <class PointTraits_, class Distance_, index_type Index_>
void DistCoverTree<PointTraits_, Distance_, Index_>::hub_query(const Point& query, Real ghost_radius, IndexVector& hub_ids) const
{
    IndexVector stack = {0};

    while (!stack.empty())
    {
        Index u = stack.back(); stack.pop_back();
        const auto& [upt, uid, uradius] = reptree[u];

        IndexVector children;
        reptree.get_children(u, children);

        auto it = ghost_map.find(uid);

        if (it != ghost_map.end())
        {
            const auto& [_, vtx] = it->second;

            if (u == vtx && distance(query, upt) <= uradius + ghost_radius)
                hub_ids.push_back(uid);
        }


        for (Index v : children)
        {
            const auto& [vpt, vid, vradius] = reptree[v];

            if (distance(query, vpt) <= vradius + ghost_radius)
                stack.push_back(v);
        }
    }
}


template <class PointTraits_, class Distance_, index_type Index_>
typename DistCoverTree<PointTraits_, Distance_, Index_>::Index
DistCoverTree<PointTraits_, Distance_, Index_>::build_replication_tree(const BallTree& repballtree)
{
    PointMap rep_point_map;
    IndexSet rep_globids_set;

    for (Index i = 0; i < repballtree.num_vertices(); ++i)
        rep_globids_set.insert(repballtree[i].id);

    IndexVector rep_globids(rep_globids_set.begin(), rep_globids_set.end());
    Index num_rep_points = rep_globids.size();

    collect_point_map(rep_globids, rep_point_map);

    auto itemizer = [&](const Ball& ball) -> PointBall { return {rep_point_map.at(ball.id), ball.id, ball.radius}; };

    repballtree.itemize_new_tree(reptree, itemizer);

    return num_rep_points;
}

template <class PointTraits_, class Distance_, index_type Index_>
void DistCoverTree<PointTraits_, Distance_, Index_>::collect_point_map(const IndexVector& globids, PointMap& point_map) const
{
        PointPairVector my_point_pairs, point_pairs;

        for (Index globid : globids)
            if (myoffset <= globid && globid < myoffset + mysize)
                my_point_pairs.emplace_back(globid, mypoints[globid-myoffset]);

        comm.allgatherv(my_point_pairs, point_pairs);

        point_map.clear();
        point_map.reserve(point_pairs.size());

        for (const auto& [globid, pt] : point_pairs)
            point_map.insert({globid, pt});
}

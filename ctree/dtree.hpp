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

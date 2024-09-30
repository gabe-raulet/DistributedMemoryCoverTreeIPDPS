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
            leaf_count += split_hub.update_tree(reptree, next_hubs, leaf_pts);
        }

        std::swap(hubs, next_hubs);

        timer.stop_timer();
        t = timer.get_max_time();
        elapsed += t;

        avg_hub_size = (totsize - leaf_count + 0.0) / hubs.size();

        if (verbose && !comm.rank())
        {
            double leaf_percent = (100.0*leaf_count)/totsize;
            fmt::print("[msg::{},elapsed={:.3f},time={:.3f}] {:.2f} percent leaves reached [iter={},levels={},vertices={},hubs={},avg_hub_size={:.3f}]\n", __func__, elapsed, t, leaf_percent, iter, reptree.num_levels(), reptree.num_vertices(), hubs.size(), avg_hub_size);
            std::cout << std::flush;
        }

        iter++;

    } while (static_cast<Real>(leaf_count) < switch_size && leaf_count < totsize);

    DistHub::free_mpi_argmax_op();
}

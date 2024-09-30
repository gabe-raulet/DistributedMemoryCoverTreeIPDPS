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
    double elapsed;
    Real avg_hub_size;
    Index leaf_count, iter, num_hubs;
    DistHubVector hubs;

    timer.start_timer();

    DistHub::init_mpi_argmax_op();

    Point root_pt = mypoints.front();
    comm.bcast(root_pt, 0);

    hubs.emplace_back(mypoints, root_pt, *this);

    timer.stop_timer();

    if (verbose && !comm.rank())
    {
        fmt::print("[time={:.3f},msg::{}] initialized root hub [farthest_point={},root_hub_radius={:.3f}]\n", timer.get_max_time(), __func__, hubs.back().cand(), hubs.back().radius());
        std::cout << std::flush;
    }

    DistHub::free_mpi_argmax_op();
}

template <class CoverTree>
Hub<CoverTree>::Hub(const PointVector& points, Index myoffset, Point repr_pt)
    : leaders({0}),
      hub_size(points.size()),
      representative(0),
      candidate(0),
      hub_parent(-1),
      hub_radius(0.),
      hub_sep(0.),
      active(true)
{
    hub_points.reserve(hub_size);

    for (Index id = 0; id < hub_size; ++id)
    {
        hub_points.emplace_back(myoffset + id, 0, distance(repr_pt, points[id]));

        if (hub_points[id].dist > hub_points[candidate].dist)
            candidate = id;
    }

    hub_radius = hub_sep = hub_points[candidate].dist;
    candidate_point = points[candidate];

    candidate += myoffset;
}

template <class CoverTree>
Hub<CoverTree>::Hub(const HubPointVector& hub_points, Index representative, Index candidate, Point candidate_point, Index parent, Real radius)
    : leaders({representative}),
      hub_points(hub_points),
      hub_size(hub_points.size()),
      representative(representative),
      candidate(candidate),
      candidate_point(candidate_point),
      hub_parent(parent),
      hub_radius(radius),
      hub_sep(radius),
      active(true) {}

template <class CoverTree>
void Hub<CoverTree>::add_new_leader(const PointVector& points)
{
    Index n = localsize();
    Index a = localoffset();

    Index new_leader = candidate;
    leaders.push_back(new_leader);

    candidate = a;
    hub_sep = 0.;

    for (Index i = 0; i < n; ++i)
    {
        auto& [id, leader, dist] = hub_points[i];
        Real new_leader_dist = distance(candidate_point, points[id-a]);

        if (new_leader_dist < dist)
        {
            dist = new_leader_dist;
            leader = new_leader;
        }

        if (dist > hub_sep)
        {
            candidate = id;
            hub_sep = dist;
        }
    }

    candidate_point = points[candidate-a]; // update candidate point
}

template <class CoverTree>
void Hub<CoverTree>::split_leaders(const PointVector& points)
{
    Index n = localsize();
    Index a = localoffset();

    for (Index leader : leaders)
    {
        Index relcand = 0;
        HubPointVector new_hub_points;

        for (Index i = 0; i < n; ++i)
        {
            if (hub_points[i].leader == leader)
            {
                new_hub_points.push_back(hub_points[i]);

                if (new_hub_points.back().dist > new_hub_points[relcand].dist)
                    relcand = new_hub_points.size()-1;
            }
        }

        Index new_hub_size = new_hub_points.size();

        Real new_radius = -1.;
        Index new_candidate = -1;
        Point new_candidate_point;

        if (new_hub_size > 0)
        {
            new_radius = new_hub_points[relcand].dist;
            new_candidate = new_hub_points[relcand].id;
            new_candidate_point = points[new_hub_points[relcand].id-a];
        }

        new_hubs.emplace_back(new_hub_points, leader, new_candidate, new_candidate_point, -1, new_radius);
    }
}

template <class CoverTree>
void Hub<CoverTree>::find_leaves(Index min_hub_size)
{
    HubVector updated_new_hubs;

    for (Index i = 0; i < new_hubs.size(); ++i)
    {
        Hub& new_hub = new_hubs[i];

        if (new_hub.size() <= min_hub_size)
            for (const HubPoint& p : new_hub.hub_points)
                leaves.push_back(p.id);
        else
            updated_new_hubs.push_back(new_hub);
    }

    std::swap(updated_new_hubs, new_hubs);
}

template <class CoverTree>
typename Hub<CoverTree>::Index
Hub<CoverTree>::update_tree(BallTree& tree, HubVector& next_hubs, std::vector<bool>& leaf_flags)
{
    assert((active));
    add_hub_vertex(tree);
    active = false;

    for (Hub& new_hub : new_hubs)
    {
        new_hub.hub_parent = hub_vertex;
        next_hubs.push_back(new_hub);
    }

    for (Index leaf : leaves)
    {
        tree.add_vertex({leaf, 0.}, hub_vertex);
        leaf_flags[leaf] = true;
    }

    return leaves.size();
}

template <class CoverTree>
typename Hub<CoverTree>::Index
Hub<CoverTree>::add_hub_vertex(BallTree& tree)
{
    assert((active));
    hub_vertex = tree.add_vertex({repr(), radius()}, parent());
    return hub_vertex;
}


template <class DistCoverTree>
void DistHub<DistCoverTree>::mpi_argmax(void *_in, void *_inout, int *len, MPI_Datatype *dtype)
{
    PointBall *in = (PointBall *)_in;
    PointBall *inout = (PointBall *)_inout;

    for (int i = 0; i < *len; ++i)
        if (in[i].radius > inout[i].radius)
            inout[i] = in[i];
}

template <class DistCoverTree>
DistHub<DistCoverTree>::DistHub(const PointVector& mypoints, Point repr_pt, const DistCoverTree& dtree) : BaseHub(mypoints, dtree.getmyoffset(), repr_pt)
{
    auto comm = dtree.getcomm();

    assert((mypoints.size() == BaseHub::size() && mypoints.size() == dtree.getmysize()));

    this->myoffset = dtree.getmyoffset();
    this->my_hub_size = dtree.getmysize();
    BaseHub::hub_size = dtree.gettotsize();

    PointBall cand_ball = BaseHub::get_cand_ball();

    comm.allreduce(cand_ball, MPI_ARGMAX);

    BaseHub::candidate_point = cand_ball.pt;
    BaseHub::candidate = cand_ball.id;
    BaseHub::hub_radius = BaseHub::hub_sep = cand_ball.radius;
}

template <class DistCoverTree>
void DistHub<DistCoverTree>::synchronize_hubs(DistHubVector& hubs, const Comm& comm)
{
    PointBallVector cand_balls(hubs.size());
    std::transform(hubs.begin(), hubs.end(), cand_balls.begin(), [](const DistHub& hub) { return hub.get_cand_ball(); });

    comm.allreduce(cand_balls, MPI_ARGMAX);

    for (Index i = 0; i < hubs.size(); ++i)
    {
        hubs[i].candidate_point = cand_balls[i].pt;
        hubs[i].candidate = cand_balls[i].id;
        hubs[i].hub_sep = cand_balls[i].radius;
    }
}

template <class DistCoverTree>
void DistHub<DistCoverTree>::synchronize_split_hubs(DistHubVector& split_hubs, Index min_hub_size, const Comm& comm)
{
    using IndexPair = std::pair<Index, Index>;
    using IndexPairVector = std::vector<IndexPair>;

    PointBallVector cand_balls;
    IndexPairVector myleaves, leaves;
    IndexVector my_new_hub_ptrs, new_hub_sizes;
    Index num_split_hubs = split_hubs.size();

    Index tot_new_hubs = 0;

    for (Index i = 0; i < num_split_hubs; ++i)
    {
        tot_new_hubs += split_hubs[i].new_hubs.size();
    }

    new_hub_sizes.reserve(tot_new_hubs);
    cand_balls.reserve(tot_new_hubs);

    for (const DistHub& split_hub : split_hubs)
        for (const BaseHub& new_hub : split_hub.new_hubs)
        {
            new_hub_sizes.push_back(new_hub.size());
            cand_balls.push_back(new_hub.get_cand_ball());
        }

    comm.allreduce(new_hub_sizes, MPI_SUM);
    comm.allreduce(cand_balls, MPI_ARGMAX);

    Index ptr = 0;

    for (Index slot = 0; slot < num_split_hubs; ++slot)
    {
        DistHub& split_hub = split_hubs[slot];

        HubVector updated_new_hubs;
        HubVector& new_hubs = split_hub.new_hubs;
        assert((split_hub.new_hub_sizes.empty()));

        for (Index i = 0; i < new_hubs.size(); ++i, ++ptr)
        {
            BaseHub& new_hub = new_hubs[i];

            if (new_hub_sizes[ptr] <= min_hub_size)
            {
                for (const HubPoint& p : new_hub.get_hub_points())
                    myleaves.emplace_back(p.id, slot);
            }
            else
            {
                updated_new_hubs.emplace_back(new_hub.get_hub_points(), new_hub.repr(), cand_balls[ptr].id, cand_balls[ptr].pt, new_hub.parent(), cand_balls[ptr].radius);
                split_hub.new_hub_sizes.push_back(new_hub_sizes[ptr]);
            }
        }

        std::swap(updated_new_hubs, split_hub.new_hubs);
    }

    comm.allgatherv(myleaves, leaves);

    for (const auto& [leaf, slot] : leaves)
    {
        split_hubs[slot].leaves.push_back(leaf);
    }
}

template <class DistCoverTree>
typename DistHub<DistCoverTree>::Index
DistHub<DistCoverTree>::update_tree(BallTree& tree, DistHubVector& next_hubs, std::vector<bool>& leaf_flags)
{
    assert((BaseHub::active));
    BaseHub::add_hub_vertex(tree);
    BaseHub::active = false;

    for (Index i = 0; BaseHub& new_hub : BaseHub::new_hubs)
    {
        new_hub.hub_parent = BaseHub::hub_vertex; // TODO: see if this is necessary
        next_hubs.emplace_back(new_hub, BaseHub::hub_vertex, new_hub_sizes[i], myoffset);
        i++;
    }

    for (Index leaf : BaseHub::leaves)
    {
        tree.add_vertex({leaf, 0.}, BaseHub::hub_vertex);
        leaf_flags[leaf] = true;
    }

    return BaseHub::leaves.size();
}

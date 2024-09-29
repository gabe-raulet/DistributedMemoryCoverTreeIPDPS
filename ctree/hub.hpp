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
Hub<CoverTree>::add_hub_vertex(BallTree& tree)
{
    assert((active));
    hub_vertex = tree.add_vertex({repr(), radius()}, parent());
    return hub_vertex;
}

template <class CoverTree>
typename Hub<CoverTree>::Index
Hub<CoverTree>::add_hub_leaves(BallTree& tree, HubVector& next_hubs)
{
    assert((active));
    active = false;

    for (Hub& new_hub : new_hubs)
    {
        new_hub.hub_parent = hub_vertex;
        next_hubs.push_back(new_hub);
    }

    for (Index leaf : leaves)
    {
        tree.add_vertex({leaf, 0.}, hub_vertex);
    }

    return leaves.size();
}

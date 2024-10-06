#ifndef HUB_H_
#define HUB_H_

#include "mpienv.h"
#include "fmt/core.h"
#include "fmt/ranges.h"
#include <assert.h>

template <class CoverTree>
class Hub
{
    public:

        using PointTraits = typename CoverTree::PointTraits;
        using Distance = typename CoverTree::Distance;
        using Index = typename CoverTree::Index;
        using Real = typename CoverTree::Real;
        using Point = typename CoverTree::Point;
        using Ball = typename CoverTree::Ball;
        using PointBall = typename CoverTree::PointBall;
        using RealVector = typename CoverTree::RealVector;
        using IndexVector = typename CoverTree::IndexVector;
        using PointVector = typename CoverTree::PointVector;
        using BallTree = typename CoverTree::BallTree;
        using PointBallTree = typename CoverTree::PointBallTree;
        using IndexSet = typename CoverTree::IndexSet;

        using BallVector = std::vector<Ball>;
        using PointBallVector = std::vector<PointBall>;

        static inline constexpr Distance distance = Distance();

        struct HubPoint
        {
            Index id, leader;
            Real dist;

            Ball ball() const { return {id, dist}; }
            PointBall point_ball(Point pt) const { return {pt, id, dist}; }
        };

        using HubVector = std::vector<Hub>;
        using HubPointVector = std::vector<HubPoint>;

        Hub(const PointVector& points) : Hub(points, 0, points.front()) {}
        Hub(const HubPointVector& hub_points, Index representative, Index candidate, Point candidate_point, Index parent, Real radius);

        Index size() const { return hub_size; }
        Index repr() const { return representative; }
        Index cand() const { return candidate; }
        Index parent() const { return hub_parent; }
        Real radius() const { return hub_radius; }
        Real sep() const { return hub_sep; }
        bool is_split(Real split_ratio) const { return sep() <= split_ratio*radius(); }
        PointBall get_cand_ball() const { return {candidate_point, candidate, hub_sep}; }
        const HubPointVector& get_hub_points() const { return hub_points; }
        const IndexVector& get_leaves() const { return leaves; }

        void add_new_leader(const PointVector& points);
        void split_leaders(const PointVector& points);
        void find_leaves(Index min_hub_size);

        Index add_hub_vertex(BallTree& tree);
        Index update_tree(BallTree& tree, HubVector& next_hubs, std::vector<bool>& leaf_flags);

    protected:

        Hub(const PointVector& points, Index myoffset, Point repr_pt);

        IndexVector leaders; /* point ids of hub leaders */
        HubPointVector hub_points; /* indices, leader pointers and leader distances of hub points */

        Index hub_size; /* number of hub points */
        Index hub_vertex; /* vertex id of hub tree vertex */
        Index hub_parent; /* vertex id of hub parent tree vertex */
        Index representative; /* point id of hub representative (first leader) */
        Index candidate; /* point id of current candidate for the next leader */
        Point candidate_point; /* current candidate point */

        Real hub_radius; /* hub radius */
        Real hub_sep; /* hub separation */

        IndexVector leaves; /* leaf point ids associated with this hub */
        HubVector new_hubs; /* child hubs */

        bool active; /* false if this hub has already contributed hub/leaf vertices to tree */

        virtual Index localsize() const { return hub_size; }
        virtual Index localoffset() const { return 0; }
};

template <class DistCoverTree>
class DistHub : public Hub<DistCoverTree>
{
    public:

        using BaseHub = Hub<DistCoverTree>;

        using typename BaseHub::PointTraits;
        using typename BaseHub::Distance;
        using typename BaseHub::Index;
        using typename BaseHub::Real;
        using typename BaseHub::Point;
        using typename BaseHub::Ball;
        using typename BaseHub::PointBall;
        using typename BaseHub::RealVector;
        using typename BaseHub::IndexVector;
        using typename BaseHub::PointVector;
        using typename BaseHub::BallTree;
        using typename BaseHub::PointBallTree;
        using typename BaseHub::HubPoint;
        using typename BaseHub::HubPointVector;
        using typename BaseHub::HubVector;
        using typename BaseHub::BallVector;
        using typename BaseHub::PointBallVector;
        using typename BaseHub::IndexSet;

        using DistHubVector = std::vector<DistHub>;
        using Comm = MPIEnv::Comm;

        inline static MPI_Op MPI_ARGMAX;

        DistHub(const PointVector& mypoints, Point repr_pt, const DistCoverTree& dtree);
        DistHub(const BaseHub& hub, Index hub_parent, Index hub_size, Index myoffset);

        static void mpi_argmax(void *_in, void *_inout, int *len, MPI_Datatype *dtype);

        static void init_mpi_argmax_op() { MPI_Op_create(&mpi_argmax, 1, &MPI_ARGMAX); }
        static void free_mpi_argmax_op() { MPI_Op_free(&MPI_ARGMAX); }

        virtual Index localsize() const override { return my_hub_size; }
        virtual Index localoffset() const override { return myoffset; }

        static void synchronize_hubs(DistHubVector& hubs, const Comm& comm);
        static void synchronize_split_hubs(DistHubVector& split_hubs, Index min_hub_size, const Comm& comm);

        Index update_tree(BallTree& tree, DistHubVector& next_hubs);

    private:

        Index my_hub_size;
        Index myoffset;
        IndexVector new_hub_sizes; // global sizes of child hubs

};

#include "hub.hpp"

#endif

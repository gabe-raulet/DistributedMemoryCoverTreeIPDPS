#ifndef DIST_COVER_TREE_H_
#define DIST_COVER_TREE_H_

#include "hub.h"
#include "ctree.h"
#include "itree.h"
#include "fmt/core.h"
#include "fmt/ranges.h"
#include "mpienv.h"
#include <assert.h>
#include <unordered_map>

template <class PointTraits_, class Distance_, index_type Index_>
class DistCoverTree
{
    public:

        using GhostTree = CoverTree<PointTraits_, Distance_, Index_>;

        using PointTraits = GhostTree::PointTraits;
        using Distance = GhostTree::Distance;
        using Index = GhostTree::Index;

        static inline constexpr Distance distance = Distance();

        using Real = GhostTree::Real;
        using Point = GhostTree::Point;

        using RealVector = GhostTree::RealVector;
        using IndexVector = GhostTree::IndexVector;
        using PointVector = GhostTree::PointVector;

        using Ball = GhostTree::Ball;
        using PointBall = GhostTree::PointBall;

        using BallTree = GhostTree::BallTree;
        using PointBallTree = GhostTree::PointBallTree;

        using IndexSet = GhostTree::IndexSet;
        using IndexMap = GhostTree::IndexMap;

        using Comm = MPIEnv::Comm;
        using GhostTreeVector = std::vector<GhostTree>;
        using GhostTreeMap = std::unordered_map<Index, GhostTree>;

        using IndexPair = std::pair<Index, Index>;
        using IndexPairMap = std::unordered_map<Index, IndexPair>;
        using IndexVectorVector = std::vector<IndexVector>;
        using IndexVectorMap = std::unordered_map<Index, IndexVector>;

        using PointMap = std::unordered_map<Index, Point>;
        using PointPair = std::pair<Index, Point>;
        using PointPairVector = std::vector<PointPair>;
        using PointPairVectorMap = std::unordered_map<Index, PointPairVector>;

        DistCoverTree(const PointVector& mypoints, const Comm& comm);

        Index getmysize() const { return mysize; }
        Index getmyoffset() const { return myoffset; }
        Index gettotsize() const { return totsize; }
        Comm getcomm() const { return comm; }
        const Point* my_point_data() const { return mypoints.data(); }

        void build(Real radius, Real split_ratio, Real switch_percent, Index min_hub_size, bool verbose = false);

        Index build_epsilon_graph(Real radius, IndexVectorVector& myneighbors) const;

        int point_owner(Index globid) const { return (std::upper_bound(offsets.begin(), offsets.end(), globid) - offsets.begin()) - 1; }
        bool owns_point(Index globid) const { return myoffset <= globid && globid < myoffset + mysize; }

    private:

        Comm comm;
        IndexVector offsets;
        Index mysize, myoffset, totsize;
        PointVector mypoints;
        IndexVector rep_leaves;

        PointBallTree reptree;
        PointMap reptree_points;
        IndexMap ghost_map, hub_sizes, hub_proc_map;
        GhostTreeMap ghost_trees;

        void point_query(const Point& query, Real radius, IndexVector& neighbors, IndexVector& ghost_hubs, Index query_hub = -1) const;

        using DistHub = DistHub<DistCoverTree>;
        using DistHubVector = typename DistHub::DistHubVector;

        using PointTriple = std::tuple<Index, Index, Point>; // hub id, point id, point
        using PointTripleVector = std::vector<PointTriple>;

        template <class RandomGen>
        double estimate_workload(const DistHub& hub, RandomGen& gen, Real split_ratio, Real sample_ratio) const
        {
            Index n = static_cast<Index>(std::floor(sample_ratio * hub.localsize()));
            if (n <= 10) return 0.;

            PointVector pts;
            pts.reserve(hub.localsize());

            for (const auto& hub_point : hub.get_hub_points())
            {
                pts.push_back(mypoints[hub_point.id-myoffset]);
            }

            std::shuffle(pts.begin(), pts.end(), gen);
            pts.resize(n);

            double t = -MPI_Wtime();
            GhostTree gt(pts);
            gt.build(split_ratio, 2, false, false);
            t += MPI_Wtime();

            return t;
        }

};

#include "dtree.hpp"

#endif

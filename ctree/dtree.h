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

        using PointMap = std::unordered_map<Index, Point>;
        using PointPair = std::pair<Index, Point>;
        using PointPairVector = std::vector<PointPair>;

        DistCoverTree(const Comm& comm) : comm(comm), mysize(0), myoffset(0), totsize(0) {}
        DistCoverTree(const PointVector& points, int root, const Comm& comm);
        DistCoverTree(const PointVector& mypoints, const Comm& comm);

        Index getmysize() const { return mysize; }
        Index getmyoffset() const { return myoffset; }
        Index gettotsize() const { return totsize; }
        Comm getcomm() const { return comm; }
        const Point* my_point_data() const { return mypoints.data(); }

        void build(Real ghost_radius, Real split_ratio, Real switch_percent, Index min_hub_size, bool verbose = false);

        Index num_rep_levels() const { return reptree.num_levels(); }
        Index num_rep_vertices() const { return reptree.num_vertices(); }

        Index build_epsilon_graph(Real radius, IndexVectorVector& myneighbors, bool verbose = false) const;
        void point_query(const Point& query, Real epsilon, IndexVector& neighbors) const;

    private:

        Comm comm;
        Index mysize, myoffset, totsize;
        PointVector mypoints;
        PointBallTree reptree;

        GhostTreeMap ghost_trees; /* (local) maps hub representative to local ghost tree (only stores hub reprs that are local) */
        IndexPairMap ghost_map; /* (global) maps hub representative to (global slot, vertex) */
        IndexMap hub_to_proc_map; /* (global) maps hub representatives to their processor owners */

        void collect_point_map(const IndexVector& globids, PointMap& point_map) const;
        Index build_replication_tree(const BallTree& repballtree);
        void hub_query(const Point& query, Real ghost_radius, IndexVector& hub_ids) const;
        void reptree_point_query(const Point& query, Real radius, IndexVector& hub_ids, IndexVector& rep_neighbors) const;
};

#include "dtree.hpp"

#endif

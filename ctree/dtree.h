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

        using CoverTree = CoverTree<PointTraits_, Distance_, Index_>;

        using PointTraits = CoverTree::PointTraits;
        using Distance = CoverTree::Distance;
        using Index = CoverTree::Index;

        static inline constexpr Distance distance = Distance();

        using Real = CoverTree::Real;
        using Point = CoverTree::Point;

        using RealVector = CoverTree::RealVector;
        using IndexVector = CoverTree::IndexVector;
        using PointVector = CoverTree::PointVector;

        using Ball = CoverTree::Ball;
        using PointBall = CoverTree::PointBall;

        using BallTree = CoverTree::BallTree;
        using PointBallTree = CoverTree::PointBallTree;

        using IndexSet = CoverTree::IndexSet;
        using IndexMap = CoverTree::IndexMap;

        using Comm = MPIEnv::Comm;
        using CoverTreeVector = std::vector<CoverTree>;
        using CoverTreeMap = std::unordered_map<Index, CoverTree>;

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

        Index build_epsilon_graph(Real radius, IndexVectorVector& myneighbors) const;

    private:

        Comm comm;
        Index mysize, myoffset, totsize;
        PointVector mypoints;
        PointBallTree reptree;

        CoverTreeMap ghost_trees; /* (local) maps hub representative to local ghost tree (only stores hub reprs that are local) */
        IndexPairMap ghost_map; /* (global) maps hub representative to (global slot, vertex) */
        IndexMap hub_to_proc_map; /* (global) maps hub representatives to their processor owners */

        void collect_point_map(const IndexVector& globids, PointMap& point_map) const;
        Index build_replication_tree(const BallTree& repballtree);
        void hub_query(const Point& query, Real ghost_radius, IndexVector& hub_ids) const;
};

#include "dtree.hpp"

#endif

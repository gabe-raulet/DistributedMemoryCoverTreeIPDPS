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

        using Comm = MPIEnv::Comm;
        using CoverTreeVector = std::vector<CoverTree>;

        DistCoverTree(const Comm& comm) : comm(comm), mysize(0), myoffset(0), totsize(0) {}
        DistCoverTree(const PointVector& points, int root, const Comm& comm);
        DistCoverTree(const PointVector& mypoints, const Comm& comm);

        Index getmysize() const { return mysize; }
        Index getmyoffset() const { return myoffset; }
        Index gettotsize() const { return totsize; }
        Comm getcomm() const { return comm; }
        const Point* my_point_data() const { return mypoints.data(); }

        void build(Real ghost_radius, Real split_ratio, Real switch_percent, Index min_hub_size, bool verbose = false);

    private:

        Comm comm;
        Index mysize, myoffset, totsize;
        PointVector mypoints;
        BallTree reptree;
};

#include "dtree.hpp"

#endif

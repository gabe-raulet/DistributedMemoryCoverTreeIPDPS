#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include "hub.h"
#include "itree.h"
#include "fmt/core.h"
#include "fmt/ranges.h"
#include <assert.h>
#include <omp.h>

template <class PointTraits_, class Distance_, index_type Index_>
class CoverTree
{
    public:

        using PointTraits = PointTraits_;
        using Distance = Distance_;
        using Index = Index_;

        using Real = typename Distance::Real;
        using Point = typename PointTraits::Point;

        using RealVector = std::vector<Real>;
        using IndexVector = std::vector<Index>;
        using PointVector = std::vector<Point>;

        CoverTree(const PointVector& points) : points(points), size(points.size()) {}

        void build(Real ghost_radius, Real split_ratio, Real switch_size, Index min_hub_size, bool level_synch, bool threaded, bool verbose = false);
        void point_query(const Point& query, Real epsilon, IndexVector& neighbors) const;
        void write_tree_file(const char *fname) const;
        bool is_correct(Real split_ratio) const;

        Index num_levels() const { return tree.num_levels(); }
        Index num_vertices() const { return tree.num_vertices(); }
        Index num_points() const { return size; }

        struct Ball { Index id; Real radius; };
        struct PointBall { Point pt; Index id; Real radius; };

        using BallTree = InsertTree<Ball, Index>;
        using PointBallTree = InsertTree<PointBall, Index>;

        static inline constexpr Distance distance = Distance();

    protected:

        Index size;
        PointVector points;
        PointBallTree tree;
        bool ghost_trees;

        void fill_point_ball_tree(const BallTree& balltree, const PointVector& points);
};

#include "ctree.hpp"

#endif

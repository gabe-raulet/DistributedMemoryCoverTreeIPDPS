#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include "hub.h"
#include "itree.h"
#include "fmt/core.h"
#include "fmt/ranges.h"
#include <assert.h>
#include <unordered_set>
#include <unordered_map>
#include <omp.h>

#ifdef LOG
#include "json.hpp"
#include <fstream>
using json = nlohmann::json;
#endif

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
        using CoverTreeVector = std::vector<CoverTree>;

        using IndexSet = std::unordered_set<Index>;
        using IndexMap = std::unordered_map<Index, Index>;
        using IndexPair = std::pair<Index, Index>;
        using IndexPairMap = std::unordered_map<Index, IndexPair>;
        using IndexVectorVector = std::vector<IndexVector>;

        CoverTree() {}
        CoverTree(const PointVector& points) : points(points) {}
        CoverTree(const PointVector& points, const IndexVector& globids) : points(points), globids(globids) {}

        void build(Real split_ratio, Index min_hub_size, bool threaded, bool verbose = false);
        void point_query(const Point& query, Real epsilon, IndexVector& neighbors) const;
        bool is_correct(Real split_ratio) const;

        Index num_levels() const { return tree.num_levels(); }
        Index num_vertices() const { return tree.num_vertices(); }
        Index num_points() const { return points.size(); }
        const Point* point_data() const { return points.data(); }

        struct Ball { Index id; Real radius; };
        struct PointBall { Point pt; Index id; Real radius; };

        using BallTree = InsertTree<Ball, Index>;
        using PointBallTree = InsertTree<PointBall, Index>;

        static inline constexpr Distance distance = Distance();

    private:

        PointVector points;
        IndexVector globids;
        PointBallTree tree;

        bool has_globids() const { return !globids.empty(); }

        void add_point(Point pt, Index globid)
        {
            points.push_back(pt);
            globids.push_back(globid);
        }

        void set_new_root(Index root)
        {
            if (has_globids())
            {
                Index where = std::find(globids.begin(), globids.end(), root) - globids.begin();

                if (where != 0)
                {
                    std::swap(globids[0], globids[where]);
                    std::swap(points[0], points[where]);
                }
            }
            else if (root != 0) std::swap(points[0], points[root]);
        }
};

template <class PointTraits_, class Distance_, index_type Index_>
class GhostTree
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
        using GhostTreeVector = std::vector<GhostTree>;

        using IndexSet = std::unordered_set<Index>;
        using IndexMap = std::unordered_map<Index, Index>;
        using IndexPair = std::pair<Index, Index>;
        using IndexPairMap = std::unordered_map<Index, IndexPair>;
        using IndexVectorVector = std::vector<IndexVector>;

        GhostTree() : size(0) {}
        GhostTree(const PointVector& points) : points(points), size(points.size()) {}
        GhostTree(const PointVector& points, const IndexVector& globids) : points(points), globids(globids), size(points.size()) {}

        void build(Real ghost_radius, Real split_ratio, Real switch_percent, Index min_hub_size, bool level_synch, bool threaded, bool verbose = false);
        void point_query(const Point& query, Real epsilon, IndexVector& neighbors) const;
        void write_tree_file(const char *fname) const;
        bool is_correct(Real split_ratio) const;

        Index num_levels() const { return tree.num_levels(); }
        Index num_vertices() const { return tree.num_vertices(); }
        Index num_points() const { return size; }
        const Point* point_data() const { return points.data(); }

        struct Ball { Index id; Real radius; };
        struct PointBall { Point pt; Index id; Real radius; };

        using BallTree = InsertTree<Ball, Index>;
        using PointBallTree = InsertTree<PointBall, Index>;

        static inline constexpr Distance distance = Distance();

    protected:

        Index size;
        PointVector points;
        BallTree tree;

        IndexVector globids;
        IndexPairMap ghost_map; /* maps hub representative to (slot, vertex) */
        GhostTreeVector ghost_trees;

        bool has_ghost_trees() const { return !ghost_trees.empty(); }
        bool has_globids() const { return !globids.empty(); }
        void hub_query(const Point& query, Real ghost_radius, IndexVector& hub_ids) const;
        void reptree_point_query(const Point& query, Real epsilon, IndexVector& neighbors) const;
        void ghost_point_query(const Point& query, Real epsilon, IndexVector& neighbors) const;

        void add_point(Point pt, Index globid)
        {
            points.push_back(pt);
            globids.push_back(globid);
            size++;
        }

        void set_new_root(Index root)
        {
            if (has_globids())
            {
                Index where = std::find(globids.begin(), globids.end(), root) - globids.begin();

                if (where != 0)
                {
                    std::swap(globids[0], globids[where]);
                    std::swap(points[0], points[where]);
                }
            }
            else if (root != 0) std::swap(points[0], points[root]);
        }


};

#include "ctree.hpp"

#endif

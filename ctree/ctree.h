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

        void build(Real split_ratio, Index min_hub_size, bool threaded, json& stats_json, bool verbose = false);
        void point_query(const Point& query, Real epsilon, IndexVector& neighbors) const;
        bool is_correct(Real split_ratio) const;

        Index num_levels() const { return tree.num_levels(); }
        Index num_vertices() const { return tree.num_vertices(); }
        Index num_points() const { return points.size(); }
        const Point* point_data() const { return points.data(); }
        const Index* globid_data() const { return globids.data(); }

        struct Ball { Index id; Real radius; };
        struct PointBall { Point pt; Index id; Real radius; };

        using BallTree = InsertTree<Ball, Index>;
        using PointBallTree = InsertTree<PointBall, Index>;

        static inline constexpr Distance distance = Distance();

        void add_point(Point pt, Index globid);
        void set_new_root(Index root);

        Index build_epsilon_graph(Real radius, IndexVectorVector& neighbors) const;

    private:

        PointVector points;
        IndexVector globids;
        PointBallTree tree;

        bool has_globids() const { return !globids.empty(); }
};

#include "ctree.hpp"

#endif

#ifndef VECTOR_POINT_TRAITS_H_
#define VECTOR_POINT_TRAITS_H_

#include <array>
#include <vector>
#include <unordered_set>
#include <concepts>
#include <algorithm>
#include <functional>
#include <fstream>
#include <assert.h>

template <int D>
concept is_pos_int = (D >= 1);

template <class Atom, int D> requires is_pos_int<D>
struct vector_point_type
{
    using Point = std::array<Atom, D>;
};

template <class Atom>
struct vector_point_type<Atom, 1>
{
    using Point = Atom;
};

template <class Atom_, int D> requires is_pos_int<D>
struct VectorPointTraits
{
    using Atom = Atom_;

    static inline constexpr int dimension = D;

    using Point = typename vector_point_type<Atom, D>::Point;
    using PointRecord = std::array<char, sizeof(int) + sizeof(Point)>;

    static inline consteval size_t dim_size() { return sizeof(D); }
    static inline consteval size_t point_size() { return sizeof(Point); }
    static inline consteval size_t record_size() { return sizeof(PointRecord); }

    struct PointHash { size_t operator()(const Point& p) const noexcept; };

    using PointVector = std::vector<Point>;
    using PointSet = std::unordered_set<Point, PointHash>;

    static void pack_point(PointRecord& record, const Point& p);
    static void unpack_point(const PointRecord& record, Point& p);

    static void read_from_file(PointVector& points, const char *fname);
    static void write_to_file(const PointVector& points, const char *fname);

    template <class RandomGen, class RandomDist>
    static void fill_random_point(Point& point, RandomGen& gen, RandomDist& dist);

    template <class RandomGen, class RandomDist>
    static void fill_random_points(PointVector& points, RandomGen& gen, RandomDist& dist);

    template <class RandomGen, class RandomDist, class Iter>
    static void fill_random_points(Iter first, Iter last, RandomGen& gen, RandomDist& dist);
};

#include "vptraits.hpp"

#endif

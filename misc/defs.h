#ifndef DIM_SIZE
#error "DIM_SIZE is undefined"
#endif

#if (FP_SIZE == 64)
using Real = double;
#else
using Real = float;
#endif

using PointTraits = VectorPointTraits<Real, DIM_SIZE>;
using Distance = L2Distance<PointTraits>;

using Point = typename PointTraits::Point;

using RealVector = std::vector<Real>;
using IndexVector = std::vector<Index>;
using PointVector = std::vector<Point>;

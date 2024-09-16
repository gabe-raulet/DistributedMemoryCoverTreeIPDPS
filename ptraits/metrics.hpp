template <class Atom_, real_type Real_, int D>
struct L2Distance<VectorPointTraits<Atom_, D>, Real_>
{
    using Atom = Atom_;
    using Real = Real_;

    using PointTraits = VectorPointTraits<Atom_, D>;
    using Point = typename PointTraits::Point;

    Real operator()(const Point& p, const Point& q) const
    {
        if constexpr (D == 1) return p < q? static_cast<Real>(q - p) : static_cast<Real>(p - q);
        else
        {
            Real sum = 0., delta;
            for (int i = 0; i < D; ++i)
            {
                delta = static_cast<Real>(p[i] - q[i]);
                sum += delta * delta;
            }
            return std::sqrt(sum);
        }
    }
};

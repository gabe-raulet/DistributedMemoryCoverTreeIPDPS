template <class Atom_, int D> requires is_pos_int<D>
size_t VectorPointTraits<Atom_, D>::PointHash::operator()(const Point& p) const noexcept
{
    if constexpr (std::same_as<Atom, Point>) return std::hash<Atom>{}(p);
    else
    {
        size_t h = 0;
        for (const Atom& e : p) h ^= std::hash<Atom>{}(e) + (h << 6) + (h >> 2);
        return h;
    }
}

template <class Atom_, int D> requires is_pos_int<D>
void VectorPointTraits<Atom_, D>::pack_point(PointRecord& record, const Point& p)
{
    char *dim_dest = record.data();
    char *pt_dest = dim_dest + dim_size();
    char *pt_src;

    if constexpr (D == 1) pt_src = (char*)&p;
    else pt_src = (char*)p.data();

    std::memcpy(dim_dest, &dimension, dim_size());
    std::memcpy(pt_dest, pt_src, point_size());
}

template <class Atom_, int D> requires is_pos_int<D>
void VectorPointTraits<Atom_, D>::unpack_point(const PointRecord& record, Point& p)
{
    int dim;

    const char *dim_src = record.data();
    const char *pt_src = record.data() + dim_size();
    char *pt_dest;

    if constexpr (D == 1) pt_dest = (char*)&p;
    else pt_dest = (char*)p.data();

    std::memcpy(&dim, dim_src, dim_size()); assert(dim == D);
    std::memcpy(pt_dest, pt_src, point_size());
}

template <class Atom_, int D> requires is_pos_int<D>
void VectorPointTraits<Atom_, D>::read_from_file(PointVector& points, const char *fname)
{
    Point p;
    std::ifstream is;
    PointRecord record;
    size_t filesize, n;
    int dim;

    is.open(fname, std::ios::binary | std::ios::in);

    is.seekg(0, is.end);
    filesize = is.tellg();
    is.seekg(0, is.beg);

    is.read((char*)&dim, dim_size()); assert(dim == D);
    is.seekg(0, is.beg);

    assert(filesize % record_size() == 0);
    n = filesize / record_size();
    points.resize(n);

    for (size_t i = 0; i < n; ++i)
    {
        is.read(record.data(), record_size());
        unpack_point(record, points[i]);
    }

    is.close();
}

template <class Atom_, int D> requires is_pos_int<D>
void VectorPointTraits<Atom_, D>::write_to_file(const PointVector& points, const char *fname)
{
    std::ofstream os;
    PointRecord record;

    os.open(fname, std::ios::binary | std::ios::out);

    for (const Point point : points)
    {
        pack_point(record, point);
        os.write(record.data(), record_size());
    }

    os.close();
}

template <class Atom_, int D> requires is_pos_int<D>
template <class RandomGen, class RandomDist>
void VectorPointTraits<Atom_, D>::fill_random_point(Point& point, RandomGen& gen, RandomDist& dist)
{
    if constexpr (D == 1) point = dist(gen);
    else std::generate(point.begin(), point.end(), [&] () { return dist(gen); });
}

template <class Atom_, int D> requires is_pos_int<D>
template <class RandomGen, class RandomDist>
void VectorPointTraits<Atom_, D>::fill_random_points(PointVector& points, RandomGen& gen, RandomDist& dist)
{
    std::for_each(points.begin(), points.end(), [&] (Point& p) { fill_random_point(p, gen, dist); });
}

template <class Atom_, int D> requires is_pos_int<D>
template <class RandomGen, class RandomDist, class Iter>
void VectorPointTraits<Atom_, D>::fill_random_points(Iter first, Iter last, RandomGen& gen, RandomDist& dist)
{
    std::for_each(first, last, [&] (Point& p) { fill_random_point(p, gen, dist); });
}

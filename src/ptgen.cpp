#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif

#include "fmt/core.h"
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <limits>
#include <unistd.h>
#include "misc.h"
#include "timer.h"
#include "vptraits.h"
#include "metrics.h"
#include "version.h"
#include <omp.h>

using Index = int64_t;

#include "defs.h"

int main(int argc, char *argv[])
{
    Index size;
    PointVector points;
    const char *fname = NULL;

    int seed = -1;
    Real var = 10.0;

    auto usage = [&] (int err)
    {
        fprintf(stderr, "Usage: %s [options] <size> <filename>\n", argv[0]);
        fprintf(stderr, "Options: -s INT   ptgen rng seed [%d]\n", seed);
        fprintf(stderr, "         -V FLOAT ptgen variance [%.2f]\n", var);
        fprintf(stderr, "         -h       help message\n");
        exit(err);
    };

    int c;
    while ((c = getopt(argc, argv, "s:V:t:b:h")) >= 0)
    {
        if      (c == 's') seed = atoi(optarg);
        else if (c == 'V') var = atof(optarg);
        else if (c == 'h') usage(0);
    }

    if (argc - optind < 2)
    {
        fmt::print(stderr, "[err::{}] missing argument(s)\n", __func__);
        usage(1);
    }

    size = read_integer<Index>(argv[optind++]);
    fname = argv[optind];

    if (seed < 0) // generate random seed if none provided
    {
        std::random_device rd;
        seed = rd();
        seed = seed < 0? -seed : seed;
    }

    std::stringstream ss;
    for (int i = 0; i < argc-1; ++i) ss << argv[i] << " "; ss << argv[argc-1];
    fmt::print("[msg::{}] [when='{}',commit='{}'] cmd: {}\n", __func__, return_current_date_and_time(), GIT_COMMIT, ss.str());

    double t;

    t = -omp_get_wtime();
    points.resize(size);
    std::default_random_engine gen(seed);
    std::normal_distribution<Real> dist{0.0, std::sqrt(var)};
    PointTraits::fill_random_points(points, gen, dist);
    t += omp_get_wtime();

    fmt::print("[time={:.3f},msg::{}] generated {} random points [fp={},dim={},var={:.2f},seed={}]\n", t, __func__, size, sizeof(Real)<<3, DIM_SIZE, var, seed);

    t = -omp_get_wtime();
    PointTraits::write_to_file(points, fname);
    t += omp_get_wtime();

    fmt::print("[time={:.3f},msg::{}] wrote {} points to file '{}'\n", t, __func__, size, fname);

    return 0;
}


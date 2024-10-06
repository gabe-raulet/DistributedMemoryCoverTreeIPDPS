#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif

#include "fmt/core.h"
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <unistd.h>
#include "misc.h"
#include "timer.h"
#include "dtree.h"
#include "vptraits.h"
#include "metrics.h"
#include "mpienv.h"
#include "version.h"
#include <omp.h>

#ifndef DIM_SIZE
#error "DIM_SIZE is undefined"
#endif

using Index = int64_t;
using Comm = MPIEnv::Comm;

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

Index mysize; /* my rank's number of points */
Index totsize; /* total number of points across all ranks */
PointVector mypoints; /* my rank's points */
const char *fname = NULL; /* input points filename */

Real radius = 0.; /* epsilon graph radius */
Real split_ratio = 0.5; /* split hubs distance ratio */
Real switch_percent = 100.; /* switch to ghost hub/tree task parallelism when we reach this perctentage of leaves */
Index min_hub_size = 10; /* hubs below this size automatically become all leaves */

bool build_graph = false;
bool verbose = false;

void parse_arguments(int argc, char *argv[], const Comm& comm);

int main_mpi(int argc, char *argv[]);
int main(int argc, char *argv[])
{
    MPIEnv::initialize(&argc, &argv);
    int err = main_mpi(argc, argv);
    MPIEnv::finalize();
    return err;
}

int main_mpi(int argc, char *argv[])
{
    auto comm = Comm::world();

    parse_arguments(argc, argv, comm);

    auto timer = comm.get_timer();
    timer.start_timer();

    {
        PointVector points;
        std::vector<int> sendcounts;

        if (!comm.rank())
        {
            PointTraits::read_from_file(points, fname);
            totsize = points.size();

            sendcounts.resize(comm.size());
            get_balanced_counts(sendcounts, (size_t)totsize);
        }

        comm.scatterv(points, sendcounts, mypoints, 0);
        mysize = mypoints.size();
    }

    timer.stop_timer();

    if (!comm.rank()) fmt::print("[msg::{},time={:.3f}] read {} points from file '{}'\n", __func__, timer.get_max_time(), totsize, fname);

    timer.start_timer();
    DistCoverTree<PointTraits, Distance, Index> dtree(mypoints, comm);
    timer.stop_timer();

    if (!comm.rank()) fmt::print("[msg::{},time={:.3f}] initialized distributed cover tree\n", __func__, timer.get_max_time());

    timer.start_timer();
    dtree.build(radius, split_ratio, switch_percent, min_hub_size, verbose);
    timer.stop_timer();

    if (!comm.rank()) fmt::print("[msg::{},time={:.3f}] constructed distributed cover tree [rep_vertices={},rep_levels={},avg_rep_nesting={:.3f}]\n", __func__, timer.get_max_time(), dtree.num_rep_vertices(), dtree.num_rep_levels(), -1.0);

    /* if (build_graph) */
    /* { */
        /* std::vector<IndexVector> mygraph; */
        /* Index num_edges; */

        /* timer.start_timer(); */
        /* num_edges = dtree.build_epsilon_graph(radius, mygraph); */
        /* timer.stop_timer(); */

        /* if (!comm.rank()) fmt::print("[msg::{},time={:.3f}] constructed epsilon graph [vertices={},edges={},avg_deg={:.3f}]\n", __func__, timer.get_max_time(), totsize, num_edges, (num_edges+0.0)/totsize); */
    /* } */

    return 0;
}

void parse_arguments(int argc, char *argv[], const Comm& comm)
{
    auto usage = [&] (int err, bool isroot)
    {
        if (isroot)
        {
            fprintf(stderr, "Usage: %s [options] <filename>\n", argv[0]);
            fprintf(stderr, "Options: -r FLOAT  graph radius [optional]\n");
            fprintf(stderr, "         -S FLOAT  hub split ratio [%.2f]\n", split_ratio);
            fprintf(stderr, "         -s FLOAT  switch percent [%.2f]\n", switch_percent);
            fprintf(stderr, "         -l INT    minimum hub size [%lu]\n", (size_t)min_hub_size);
            fprintf(stderr, "         -v        verbose\n");
            fprintf(stderr, "         -h        help message\n");
        }

        MPIEnv::exit(err);
    };

    int c;
    while ((c = getopt(argc, argv, "r:S:s:l:vh")) >= 0)
    {
        if      (c == 'r') radius = atof(optarg);
        else if (c == 'S') split_ratio = atof(optarg);
        else if (c == 's') switch_percent = atof(optarg);
        else if (c == 'l') min_hub_size = atoi(optarg);
        else if (c == 'v') verbose = true;
        else if (c == 'h') usage(0, !comm.rank());
    }

    if (argc - optind < 1)
    {
        if (!comm.rank()) fmt::print(stderr, "[err::{}] missing argument(s)\n", __func__);
        usage(1, !comm.rank());
    }

    fname = argv[optind];
    build_graph = (radius > 0);

    if (!comm.rank())
    {
        std::stringstream ss;
        for (int i = 0; i < argc-1; ++i) ss << argv[i] << " "; ss << argv[argc-1];
        fmt::print("[msg::{},mpi_num_ranks={},commit={},when='{}'] cmd: {}\n", __func__, comm.size(), GIT_COMMIT, return_current_date_and_time(), ss.str());
        fmt::print("[msg::{}] point parameters: [file='{}',dim={},fp={}]\n", __func__, fname, DIM_SIZE, sizeof(Real)<<3);
        fmt::print("[msg::{}] ctree parameters: [split_ratio={:.2f},switch_percent={:.2f},min_hub_size={},verbose={}]\n", __func__, split_ratio, switch_percent, min_hub_size, verbose);
        if (build_graph) fmt::print("[msg::{}] graph parameters: [radius={:.3f}]\n", __func__, radius);
    }
}

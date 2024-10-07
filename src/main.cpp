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
#include "ctree.h"
#include "vptraits.h"
#include "metrics.h"
#include "version.h"
#include <omp.h>

#ifndef DIM_SIZE
#error "DIM_SIZE is undefined"
#endif

using Index = int64_t;

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

Index size; /* total number of points */
PointVector points; /* points */
const char *fname = NULL; /* input points filename */
const char *stats_fname = "output.json"; /* output json filename with statistics */

Real radius = 0.; /* epsilon graph radius (radius <= 0 means that we don't build epsilon graph) */
Real split_ratio = 0.5; /* split hubs distance ratio */
Real switch_percent = 100.0; /* switch to ghost hub/tree task parallelism when we reach this percentage of leaves */
Index min_hub_size = 10; /* hubs below this size automatically become all leaves */

bool verify_tree = false; /* verify correctness of constructed cover tree */
bool verify_graph = false; /* verify correctness of constructed epsilon neighbor graph (slow because it builds brute force graph to compare) */
bool build_graph = false;
bool verbose = false;

void parse_arguments(int argc, char *argv[]);
bool graph_is_correct(const PointVector& points, Real radius, const std::vector<IndexVector>& graph);

int main(int argc, char *argv[])
{
    double t;

    parse_arguments(argc, argv);

    t = -omp_get_wtime();
    PointTraits::read_from_file(points, fname); size = points.size();
    t += omp_get_wtime();

    fmt::print("[msg::{},time={:.3f}] read {} points from file '{}'\n", __func__, t, size, fname);

    int nthreads;

    #pragma omp parallel
    nthreads = omp_get_num_threads();

    json stats_json;
    json parameter_info, points_info, build_info;

    parameter_info["radius"] = radius;
    parameter_info["split_ratio"] = split_ratio;
    parameter_info["switch_percent"] = switch_percent;
    parameter_info["min_hub_size"] = min_hub_size;
    parameter_info["num_threads"] = nthreads;

    points_info["filename"] = fname;
    points_info["size"] = size;
    points_info["dimension"] = DIM_SIZE;
    points_info["float_size"] = sizeof(Real) << 3;

    stats_json["points_info"] = points_info;
    stats_json["parameter_info"] = parameter_info;

    t = -omp_get_wtime();
    CoverTree<PointTraits, Distance, Index> ctree(points);
    ctree.build(split_ratio, min_hub_size, true, build_info, verbose);
    t += omp_get_wtime();

    build_info["total_time"] = t;
    stats_json["build_info"] = build_info;

    fmt::print("[msg::{},time={:.3f}] constructed cover tree [vertices={},levels={},avg_nesting={:.3f}]\n", __func__, t, ctree.num_vertices(), ctree.num_levels(), (ctree.num_vertices()+0.0)/size);

    if (verify_tree)
    {
        t = -omp_get_wtime();
        bool passed = ctree.is_correct(split_ratio);
        t += omp_get_wtime();

        fmt::print("[msg::{},time={:.3f}] cover tree {} verification\n", __func__, t, passed? "PASSED" : "FAILED");
    }

    if (build_graph)
    {
        std::vector<IndexVector> graph;
        Index num_edges;

        t = -omp_get_wtime();
        num_edges = ctree.build_epsilon_graph(radius, graph);
        t += omp_get_wtime();

        fmt::print("[msg::{},time={:.3f}] constructed epsilon graph [vertices={},edges={},avg_deg={:.3f}]\n", __func__, t, size, num_edges, (num_edges+0.0)/size);

        if (verify_graph)
        {
            t = -omp_get_wtime();
            bool correct = graph_is_correct(points, radius, graph);
            t += omp_get_wtime();

            fmt::print("[msg::{},time={:.3f}] epsilon graph {} verification\n",  __func__, t, correct? "PASSED" : "FAILED");
        }
    }

    std::ofstream f(stats_fname);
    f << std::setw(4) << stats_json << std::endl;
    f.close();

    return 0;
}

void parse_arguments(int argc, char *argv[])
{
    auto usage = [&] (int err)
    {
        fprintf(stderr, "Usage: %s [options] <filename>\n", argv[0]);
        fprintf(stderr, "Options: -r FLOAT  graph radius [optional]\n");
        fprintf(stderr, "         -S FLOAT  hub split ratio [%.2f]\n", split_ratio);
        fprintf(stderr, "         -s FLOAT  switch percent [%.2f]\n", switch_percent);
        fprintf(stderr, "         -l INT    minimum hub size [%lu]\n", (size_t)min_hub_size);
        fprintf(stderr, "         -o FILE   output json ['%s']\n", stats_fname);
        fprintf(stderr, "         -T        verify tree correctness\n");
        fprintf(stderr, "         -G        verify graph correctness [assumes -r]\n");
        fprintf(stderr, "         -v        verbose\n");
        fprintf(stderr, "         -h        help message\n");
        exit(err);
    };

    int c;
    while ((c = getopt(argc, argv, "r:S:s:e:l:o:TGvh")) >= 0)
    {
        if      (c == 'r') radius = atof(optarg);
        else if (c == 'S') split_ratio = atof(optarg);
        else if (c == 's') switch_percent = atof(optarg);
        else if (c == 'l') min_hub_size = atoi(optarg);
        else if (c == 'o') stats_fname = optarg;
        else if (c == 'T') verify_tree = true;
        else if (c == 'G') verify_graph = true;
        else if (c == 'v') verbose = true;
        else if (c == 'h') usage(0);
    }

    if (argc - optind < 1)
    {
        fmt::print(stderr, "[err::{}] missing argument(s)\n", __func__);
        usage(1);
    }

    fname = argv[optind];
    build_graph = (radius > 0);

    int nthreads;

    #pragma omp parallel
    nthreads = omp_get_num_threads();

    std::stringstream ss;
    for (int i = 0; i < argc-1; ++i) ss << argv[i] << " "; ss << argv[argc-1];
    fmt::print("[msg::{}] cmd: {} [omp_num_threads={},commit={},when='{}']\n", __func__, ss.str(), nthreads, GIT_COMMIT, return_current_date_and_time());
    fmt::print("[msg::{}] point parameters: [file='{}',dim={},fp={}]\n", __func__, fname, DIM_SIZE, sizeof(Real)<<3);
    fmt::print("[msg::{}] ctree parameters: [split_ratio={:.2f},switch_percent={:.2f},min_hub_size={},verify_tree={},verbose={}]\n", __func__, split_ratio, switch_percent, min_hub_size, verify_tree, verbose);
    if (build_graph) fmt::print("[msg::{}] graph parameters: [radius={:.3f},verify_graph={}]\n", __func__, radius, verify_graph);
}

bool graph_is_correct(const PointVector& points, Real radius, const std::vector<IndexVector>& graph)
{
    auto distance = Distance();
    bool correct = true;

    Index size = points.size();

    #pragma omp parallel for
    for (Index i = 0; i < size; ++i)
    {
        if (!correct) continue;

        IndexVector neighbors; neighbors.reserve(graph[i].size());

        for (Index j = 0; j < size; ++j)
            if (distance(points[i], points[j]) <= radius)
                neighbors.push_back(j);

        if (neighbors.size() != graph[i].size() || !std::is_permutation(neighbors.begin(), neighbors.end(), graph[i].begin()))
        {
            #pragma omp atomic write
            correct = false;
        }
    }

    return correct;
}

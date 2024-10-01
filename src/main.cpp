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

Real radius = 0.; /* epsilon graph radius (radius <= 0 means that we don't build epsilon graph) */
Real split_ratio = 0.5; /* split hubs distance ratio */
Real switch_percent = 100.0; /* switch to ghost hub/tree task parallelism when we reach this percentage of leaves */
Index min_hub_size = 10; /* hubs below this size automatically become all leaves */

bool level_synch = true; /* level synchronous construction */
bool verify_tree = false; /* verify correctness of constructed cover tree */
bool verify_graph = false; /* verify correctness of constructed epsilon neighbor graph (slow because it builds brute force graph to compare) */
bool verify_graph_plus = false; /* count number of wrong edges found */
bool build_graph = false;
bool verbose = false;
int nthreads = 1;

void parse_arguments(int argc, char *argv[]);
bool graph_is_correct(const PointVector& points, Real radius, const std::vector<IndexVector>& graph);
std::pair<Index, Index> count_missing(const PointVector& points, Real radius, const std::vector<IndexVector>& graph);

int main(int argc, char *argv[])
{
    double t;

    parse_arguments(argc, argv);

    t = -omp_get_wtime();
    PointTraits::read_from_file(points, fname); size = points.size();
    t += omp_get_wtime();

    fmt::print("[msg::{},time={:.3f}] read {} points from file '{}'\n", __func__, t, size, fname);


    t = -omp_get_wtime();
    CoverTree<PointTraits, Distance, Index> ctree(points);
    ctree.build(split_ratio, min_hub_size, true, verbose);
    t += omp_get_wtime();

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

        if (verify_graph_plus)
        {
            t = -omp_get_wtime();
            const auto& [missing, total] = count_missing(points, radius, graph);
            t += omp_get_wtime();

            fmt::print("[msg::{},time={:.3f}] {:.3f} percent of edges missing [num_missing={},num_total={}]\n",  __func__, t, (100.*missing)/total, missing, total);
        }
        else if (verify_graph)
        {
            t = -omp_get_wtime();
            bool correct = graph_is_correct(points, radius, graph);
            t += omp_get_wtime();

            fmt::print("[msg::{},time={:.3f}] epsilon graph {} verification\n",  __func__, t, correct? "PASSED" : "FAILED");
        }
    }

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
        fprintf(stderr, "         -t INT    number of threads [%d]\n", nthreads);
        fprintf(stderr, "         -A        asynchronous tree construction\n");
        fprintf(stderr, "         -T        verify tree correctness\n");
        fprintf(stderr, "         -G        verify graph correctness [assumes -r]\n");
        fprintf(stderr, "         -C        count number of wrong edges [assumes -r]\n");
        fprintf(stderr, "         -v        verbose\n");
        fprintf(stderr, "         -h        help message\n");
        exit(err);
    };

    int c;
    while ((c = getopt(argc, argv, "r:S:s:Ct:e:l:TAGvh")) >= 0)
    {
        if      (c == 'r') radius = atof(optarg);
        else if (c == 'S') split_ratio = atof(optarg);
        else if (c == 's') switch_percent = atof(optarg);
        else if (c == 'l') min_hub_size = atoi(optarg);
        else if (c == 't') nthreads = atoi(optarg);
        else if (c == 'A') level_synch = false;
        else if (c == 'T') verify_tree = true;
        else if (c == 'G') verify_graph = true;
        else if (c == 'C') verify_graph_plus = true;
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

    omp_set_num_threads(nthreads);

    #pragma omp parallel
    nthreads = omp_get_num_threads();

    std::stringstream ss;
    for (int i = 0; i < argc-1; ++i) ss << argv[i] << " "; ss << argv[argc-1];
    fmt::print("[msg::{}] cmd: {} [omp_num_threads={},commit={},when='{}']\n", __func__, ss.str(), nthreads, GIT_COMMIT, return_current_date_and_time());
    fmt::print("[msg::{}] point parameters: [file='{}',dim={},fp={}]\n", __func__, fname, DIM_SIZE, sizeof(Real)<<3);
    fmt::print("[msg::{}] ctree parameters: [split_ratio={:.2f},switch_percent={:.2f},min_hub_size={},level_synch={},verify_tree={},verbose={}]\n", __func__, split_ratio, switch_percent, min_hub_size, level_synch, verify_tree, verbose);
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

std::pair<Index, Index> count_missing(const PointVector& points, Real radius, const std::vector<IndexVector>& graph)
{
    using IndexSet = std::unordered_set<Index>;

    auto distance = Distance();

    Index size = points.size();
    Index missing = 0;
    Index total = 0;

    #pragma omp parallel for reduction(+:missing,total)
    for (Index i = 0; i < size; ++i)
    {
        IndexVector neighbors; neighbors.reserve(graph[i].size());

        for (Index j = 0; j < size; ++j)
            if (distance(points[i], points[j]) <= radius)
                neighbors.push_back(j);

        total += neighbors.size();

        if (neighbors.size() != graph[i].size() || !std::is_permutation(neighbors.begin(), neighbors.end(), graph[i].begin()))
        {
            IndexSet incorrect(graph[i].begin(), graph[i].end());

            for (Index v : neighbors)
                if (!incorrect.contains(v))
                    missing++;
        }
    }

    return {missing, total};
}

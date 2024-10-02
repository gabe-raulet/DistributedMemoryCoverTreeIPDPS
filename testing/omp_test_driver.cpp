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

#include "json.hpp"
#include <fstream>

#ifndef DIM_SIZE
#error "DIM_SIZE is undefined"
#endif

using Index = int64_t;
using Real = float;
using PointTraits = VectorPointTraits<Real, DIM_SIZE>;
using Distance = L2Distance<PointTraits>;
using Point = typename PointTraits::Point;
using RealVector = std::vector<Real>;
using IndexVector = std::vector<Index>;
using PointVector = std::vector<Point>;

Index size;
PointVector points;

void driver(const json& test_case, json& output, const std::vector<int>& thread_counts, int trials);

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <input.json> <output.json>\n", argv[0]);
        return 1;
    }

    std::ifstream inf;
    std::ofstream outf;

    json input, output;

    inf.open(argv[1]);
    input = json::parse(inf);
    inf.close();

    double t;
    std::string fname = input["points_filename"];
    const std::vector<int> thread_counts = input["thread_counts"];
    const int trials = input["trials"];

    t = -omp_get_wtime();
    PointTraits::read_from_file(points, fname.c_str());
    size = points.size();
    t += omp_get_wtime();

    output["read_file_time"] = t;
    output["dimension"] = DIM_SIZE;
    output["float_size"] = sizeof(Real) << 3;
    output["num_points"] = size;

    std::vector<json> cases;

    for (const json& test_case : input["cases"])
    {
        cases.emplace_back();
        cases.back()["tag"] = test_case["tag"];

        driver(test_case, cases.back(), thread_counts, trials);
    }

    output["cases"] = cases;

    outf.open(argv[2]);
    outf << std::setw(4) << output << std::endl;
    outf.close();

    return 0;
}

void driver(const json& test_case, json& output, const std::vector<int>& thread_counts, int trials)
{
    Real radius = test_case["radius"];
    Real split_ratio = test_case["split_ratio"];
    Index min_hub_size = test_case["min_hub_size"];

    std::vector<json> thread_trials;

    for (int nthreads : thread_counts)
    {
        omp_set_num_threads(nthreads);

        double t;
        std::vector<double> cover_tree_times, epsilon_graph_times;
        std::vector<Index> edge_counts;

        for (int i = 0; i < trials; ++i)
        {
            t = -omp_get_wtime();
            CoverTree<PointTraits, Distance, Index> ctree(points);
            ctree.build(split_ratio, min_hub_size, true, false);
            t += omp_get_wtime();

            cover_tree_times.push_back(t);

            std::vector<IndexVector> graph;
            Index num_edges;

            t = -omp_get_wtime();
            num_edges = ctree.build_epsilon_graph(radius, graph);
            t += omp_get_wtime();

            epsilon_graph_times.push_back(t);
            edge_counts.push_back(num_edges);
        }

        thread_trials.emplace_back();
        thread_trials.back()["cover_tree_times"] = cover_tree_times;
        thread_trials.back()["epsilon_graph_times"] = epsilon_graph_times;
        thread_trials.back()["edge_counts"] = edge_counts;
    }

    output["thread_trials"] = thread_trials;
}

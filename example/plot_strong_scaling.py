
import json
import matplotlib.pyplot as plt
import numpy as np
from collections import namedtuple
import json
import matplotlib.pyplot as plt
import numpy as np
from collections import namedtuple

BuildInfo = namedtuple("BuildInfo", ["total_time", "num_hubs", "num_leaves", "num_levels", "num_vertices", "iter_times"])
GraphInfo = namedtuple("GraphInfo", ["total_time", "num_edges", "num_vertices", "dist_comps", "thread_edges", "thread_queries"])
ParameterInfo = namedtuple("ParameterInfo", ["min_hub_size", "num_threads", "radius", "split_ratio", "switch_percent"])
PointsInfo = namedtuple("PointsInfo", ["dimension", "float_size", "filename", "size"])

def parse_build_info(build_info):
    total_time = build_info["total_time"]
    iterations = build_info["iterations"]
    num_hubs = np.array([it["num_hubs"] for it in iterations])
    num_leaves = np.array([it["num_leaves"] for it in iterations])
    num_levels = np.array([it["num_levels"] for it in iterations])
    num_vertices = np.array([it["num_vertices"] for it in iterations])
    iter_times = np.array([it["time"] for it in iterations])
    return BuildInfo(total_time=total_time, num_hubs=num_hubs, num_leaves=num_leaves, num_levels=num_levels, num_vertices=num_vertices, iter_times=iter_times)

def parse_graph_info(graph_info):
    total_time = graph_info["total_time"]
    graph_threads = sorted(graph_info["graph_threads"], key=lambda o: o["thread_id"])
    num_edges = graph_info["num_edges"]
    num_vertices = graph_info["num_vertices"]
    dist_comps = np.array([t["thread_dist_comps"] for t in graph_threads])
    thread_edges = np.array([t["thread_edges"] for t in graph_threads])
    thread_queries = np.array([t["thread_queries"] for t in graph_threads])
    thread_times = np.array([t["thread_time"] for t in graph_threads])
    return GraphInfo(total_time=total_time, num_edges=num_edges, num_vertices=num_vertices, dist_comps=dist_comps, thread_edges=thread_edges, thread_queries=thread_queries)

def parse_parameter_info(parameter_info):
    return ParameterInfo(min_hub_size=parameter_info["min_hub_size"], num_threads=parameter_info["num_threads"], radius=parameter_info["radius"], split_ratio=parameter_info["split_ratio"], switch_percent=parameter_info["switch_percent"])

def parse_points_info(points_info):
    return PointsInfo(dimension=points_info["dimension"], float_size=points_info["float_size"], filename=points_info["filename"], size=points_info["size"])

def get_points_info_string(pi):
    return f"{pi.dimension}d.{pi.float_size}fp.{pi.size}size"

def get_params_info_string(pi):
    return f"{pi.min_hub_size}l.{pi.num_threads}t.{pi.radius:.2f}r.{pi.split_ratio:.2f}"

experiments = json.load(open("results.json", "r"))["experiments"]

parameter_infos = [parse_parameter_info(exp["parameter_info"]) for exp in experiments]
build_infos = [parse_build_info(exp["build_info"]) for exp in experiments]
points_infos = [parse_points_info(exp["points_info"]) for exp in experiments]
graph_infos = [parse_graph_info(exp["graph_info"]) for exp in experiments]

thread_counts = np.array([b.num_threads for b in parameter_infos])
build_times = np.array([b.total_time for b in build_infos])
graph_times = np.array([b.total_time for b in graph_infos])
combine_times = build_times + graph_times

order = np.argsort(thread_counts)
thread_counts = thread_counts[order]
build_times = build_times[order]
graph_times = graph_times[order]
combine_times = combine_times[order]

#  plt.plot(thread_counts, build_times, ':x', markersize=1.8, linewidth=0.8, color='red', label="build time")
#  plt.plot(thread_counts, graph_times, '-o', markersize=1.4, linewidth=0.8, color='blue', label="graph time")
#  plt.plot(thread_counts, combine_times, '--^', markersize=1.4, linewidth=0.8, color='orange', label="total time")
#  plt.xlabel("# of threads")
#  plt.ylabel("seconds")
#  plt.legend(loc="upper right")
#  plt.show()


build_speedups = build_times[0]/build_times
graph_speedups = graph_times[0]/graph_times
combine_speedups = combine_times[0]/combine_times

build_efficiencies = build_speedups/thread_counts
graph_efficiencies = graph_speedups/thread_counts
combine_efficiencies = combine_speedups/thread_counts

plt.plot(thread_counts, build_speedups, ':x', markersize=1.8, linewidth=0.8, color='red', label="build speedup")
plt.plot(thread_counts, graph_speedups, ':o', markersize=1.8, linewidth=0.8, color='blue', label="graph speedup")
plt.plot(thread_counts, combine_speedups, ':^', markersize=1.8, linewidth=0.8, color='green', label="total speedup")

#  plt.plot(thread_counts, build_efficiencies, '--x', markersize=1.8, linewidth=0.8, color='red', label="build efficiency")
#  plt.plot(thread_counts, graph_efficiencies, '--o', markersize=1.8, linewidth=0.8, color='blue', label="graph efficiency")
#  plt.plot(thread_counts, combine_efficiencies, '--^', markersize=1.8, linewidth=0.8, color='green', label="total efficiency")

plt.xlabel("# of threads")
plt.ylabel("speedup")
plt.legend(loc="upper left")
plt.show()


import json
import numpy as np
import matplotlib.pyplot as plt

input_json = json.load(open("v3/input.json", "r"))
output_json = json.load(open("v3/output.json", "r"))

input_cases = input_json["cases"]
output_cases = output_json["cases"]
dimension = output_json["dimension"]
float_size = output_json["float_size"]
num_points = output_json["num_points"]
read_file_time = output_json["read_file_time"]

def process_json(input_json, output_json):

    input_cases = input_json["cases"]
    output_cases = output_json["cases"]
    thread_counts = input_json["thread_counts"]

    x = np.arange(len(thread_counts))
    width = 0.20
    multiplier = 0
    fig, ax = plt.subplots()
    runtimes = dict()

    for input_case, output_case in zip(input_cases, output_cases):
        trials = output_case["thread_trials"]
        split_ratio = input_case["split_ratio"]
        num_trials = len(trials)
        edge_counts = set()
        total_times = []

        for trial in reversed(trials):
            total_times.append(np.array(trial["cover_tree_times"]) + np.array(trial["epsilon_graph_times"]))

        total_times = np.vstack(total_times)
        total_times_mean = total_times.mean(axis=1).round(3)
        runtimes[f"Split ratio = {split_ratio:.2f}"] = total_times_mean

    for attribute, measurement in runtimes.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute)
        #  ax.bar_label(rects, padding=3)
        multiplier += 1

    ax.set_ylabel("Time (seconds)")
    ax.set_xlabel("Number of threads")
    ax.set_title("OpenMP Scaling")
    ax.set_xticks(x + width, [str(tc) for tc in reversed(thread_counts)])
    ax.legend(loc="upper left")
    ax.set_ylim(0, 20)
    plt.show()

process_json(input_json, output_json)

def process_case(case, title="OpenMP Scaling"):

    trials = case["thread_trials"]
    num_trials = len(trials)
    edge_counts = set()
    cover_tree_times, epsilon_graph_times, thread_counts = [], [], []

    for trial in reversed(trials):
        cover_tree_times.append(np.array(trial["cover_tree_times"]))
        epsilon_graph_times.append(np.array(trial["epsilon_graph_times"]))
        thread_counts.append(trial["num_threads"])
        edge_counts.update(trial["edge_counts"])

    assert len(edge_counts) == 1

    cover_tree_times = np.vstack(cover_tree_times)
    epsilon_graph_times = np.vstack(epsilon_graph_times)
    total_times = cover_tree_times + epsilon_graph_times

    cover_tree_mean_times = cover_tree_times.mean(axis=1).round(2)
    epsilon_graph_mean_times = epsilon_graph_times.mean(axis=1).round(2)
    total_mean_times = total_times.mean(axis=1).round(2)

    cover_tree_std_times = cover_tree_times.std(axis=1).round(2)
    epsilon_graph_std_times = epsilon_graph_times.std(axis=1).round(2)
    total_std_times = total_times.std(axis=1).round(2)

    x = np.arange(len(thread_counts))
    width = 0.20
    multiplier = 0
    fig, ax = plt.subplots()

    runtimes = {
            "Cover Tree Construction" : cover_tree_mean_times,
            "Epsilon Graph Construction" : epsilon_graph_mean_times
    }

    for attribute, measurement in runtimes.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute)
        ax.bar_label(rects, padding=3)
        multiplier += 1

    ax.set_ylabel("Time (seconds)")
    ax.set_xlabel("Number of threads")
    ax.set_title(title)
    ax.set_xticks(x + width, [str(tc) for tc in thread_counts])
    ax.legend(loc="upper left")
    ax.set_ylim(0, 20)
    plt.show()

for i in range(len(output_cases)):
    process_case(output_cases[i], f"OpenMP scaling (Min hub size={input_cases[i]['min_hub_size']}, Split ratio={input_cases[i]['split_ratio']:.3f})")
    plt.show()

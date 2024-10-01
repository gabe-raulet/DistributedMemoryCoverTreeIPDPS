import sys
import json
from pathlib import Path
import subprocess as sp

test_json = json.load(open("example.json", "r"))
objects = test_json["objects"]
counts = test_json["counts"]
tests = test_json["tests"]

procs = []
for obj in objects:
    command = ["make", "-j3", f"D={obj['FP_SIZE']}", f"D={obj['DIM_SIZE']}", f"EXE={obj['EXE']}"]
    procs.append(sp.Popen(command))

for proc in procs:
    proc.wait()

procs = []
for count in counts:
    for obj in objects:
        command = [f"./ptgen.{obj['EXE']}", f"{count}", f"points.{obj['EXE']}.{count}"]
        procs.append(sp.Popen(command))

for proc in procs:
    proc.wait()

header=f"""#!/bin/bash

#SBATCH --nodes=4
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -o output.distributed.txt
#SBATCH --mail-user=ghraulet@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 30

"""

with open("job.sh", "w") as f:
    f.write(header)
    for count in counts:
        for obj in objects:
            points_filename = f"points.{obj['EXE']}.{count}"
            for test in tests:
                N, tpn, r, S, s, l, v = test["nodes"], test["tasks_per_node"], test["radius"], test["split_ratio"], test["switch_percent"], test["min_hub_size"], test["verbose"]
                n = N*tpn
                c = 2 * (128 // tpn)
                command = f"srun -N {N} -n {n} -c {c} --cpu_bind=cores ./main_mpi.{obj['EXE']} -v -s {s} -S {S} -l {l} -r {r} {points_filename} > output.n{n}.N{N}.{obj['EXE']}.{count}.txt"
                f.write(command + "\n")

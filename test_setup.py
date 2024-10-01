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
                command = f"srun -N {N} -n {n} -c {c} --cpu_bind=cores ./main_mpi.{obj['EXE']} -v -s {s} -S {S} -l {l} -r {r} > output.n{n}.N{N}.{obj['EXE']}.{count}.txt"
                f.write(command + "\n")



#  {"nodes" : 1, "tasks_per_node" : 8, "radius" : 0.1, "split_ratio" : 0.7, "switch_percent" : 10, "min_hub_size" : 10, "verbose" : true},
#  {"nodes" : 2, "tasks_per_node" : 8, "radius" : 0.1, "split_ratio" : 0.7, "switch_percent" : 10, "min_hub_size" : 10, "verbose" : true},
#  {"nodes" : 4, "tasks_per_node" : 8, "radius" : 0.1, "split_ratio" : 0.7, "switch_percent" : 10, "min_hub_size" : 10, "verbose" : true}

#!/bin/bash

#SBATCH --nodes=4
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -o output.distributed.txt
#SBATCH --mail-user=ghraulet@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 20

#  srun -N 1 -n 16 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.64M > output.n16.N1.s10.D3.FP64.64M.txt
#  srun -N 2 -n 32 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.64M > output.n32.N2.s10.D3.FP64.64M.txt
#  srun -N 3 -n 48 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.64M > output.n48.N3.s10.D3.FP64.64M.txt
#  srun -N 4 -n 64 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.64M > output.n64.N4.s10.D3.FP64.64M.txt


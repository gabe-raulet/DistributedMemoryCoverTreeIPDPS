#!/bin/bash

#SBATCH --nodes=1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -o output.threaded.txt
#SBATCH --mail-user=ghraulet@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 10

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

numactl -i 0-7 ./main.D3.FP64 -v -t128 -r 0.05 points.D3.FP64.1M > output.t128.D3.FP64.1M.txt
numactl -i 0-7 ./main.D3.FP64 -v -t64 -r 0.05 points.D3.FP64.1M > output.t64.D3.FP64.1M.txt
numactl -i 0-7 ./main.D3.FP64 -v -t32 -r 0.05 points.D3.FP64.1M > output.t32.D3.FP64.1M.txt
numactl -i 0-7 ./main.D3.FP64 -v -t16 -r 0.05 points.D3.FP64.1M > output.t16.D3.FP64.1M.txt
numactl -i 0-7 ./main.D3.FP64 -v -t8 -r 0.05 points.D3.FP64.1M > output.t8.D3.FP64.1M.txt
numactl -i 0-7 ./main.D3.FP64 -v -t4 -r 0.05 points.D3.FP64.1M > output.t4.D3.FP64.1M.txt
numactl -i 0-7 ./main.D3.FP64 -v -t2 -r 0.05 points.D3.FP64.1M > output.t2.D3.FP64.1M.txt
numactl -i 0-7 ./main.D3.FP64 -v -t1 -r 0.05 points.D3.FP64.1M > output.t1.D3.FP64.1M.txt

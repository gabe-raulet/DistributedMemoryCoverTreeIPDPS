#!/bin/bash

#SBATCH --nodes=4
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -o output.distributed.txt
#SBATCH --mail-user=ghraulet@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 20

srun -N 1 -n 16 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.64M > output.n16.N1.s10.D3.FP64.64M.txt
srun -N 2 -n 32 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.64M > output.n32.N2.s10.D3.FP64.64M.txt
srun -N 3 -n 48 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.64M > output.n48.N3.s10.D3.FP64.64M.txt
srun -N 4 -n 64 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.64M > output.n64.N4.s10.D3.FP64.64M.txt

#srun -N 1 -n 16 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 20 -r 0.05 points.D3.FP64.4M > output.n16.N1.s20.D3.FP64.4M.txt
#srun -N 2 -n 32 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 20 -r 0.05 points.D3.FP64.4M > output.n32.N2.s20.D3.FP64.4M.txt
#srun -N 3 -n 48 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 20 -r 0.05 points.D3.FP64.4M > output.n48.N3.s20.D3.FP64.4M.txt
#srun -N 4 -n 64 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 20 -r 0.05 points.D3.FP64.4M > output.n64.N4.s20.D3.FP64.4M.txt

#srun -n 128 -c 2 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.1M > output.n128.N1.s10.D3.FP64.1M.txt
#srun -n 64 -c 4 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.1M > output.n64.N1.s10.D3.FP64.1M.txt
#srun -n 32 -c 8 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.1M > output.n32.N1.s10.D3.FP64.1M.txt
#srun -n 16 -c 16 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.1M > output.n16.N1.s10.D3.FP64.1M.txt
#srun -n 8 -c 32 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.1M > output.n8.N1.s10.D3.FP64.1M.txt
#srun -n 4 -c 64 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.1M > output.n4.N1.s10.D3.FP64.1M.txt
#srun -n 2 -c 128 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.1M > output.n2.N1.s10.D3.FP64.1M.txt
#srun -n 1 -c 256 --cpu_bind=cores ./main_mpi.D3.FP64 -v -s 10 -r 0.05 points.D3.FP64.1M > output.n1.N1.s10.D3.FP64.1M.txt

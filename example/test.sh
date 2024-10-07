#!/bin/bash

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

export FILENAME=points
export RADIUS=2.0
export SPLIT_RATIO=0.8

for THREAD_COUNT in 1 2 4 8 16 32 64 128
do
    export OMP_NUM_THREADS=$THREAD_COUNT
    export OUTPUT=log.${THREAD_COUNT}threads.${RADIUS}radius.${SPLIT_RATIO}split_ratio.json

    numactl -i0-7 ./main -r $RADIUS -S $SPLIT_RATIO -o $OUTPUT $FILENAME
done

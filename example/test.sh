#!/bin/bash

# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads

# export FILENAME=../points
export FILENAME=sift_base_1000000.fvecs
export RADIUS=0
export SPLIT_RATIO=0.65

counter=1

for THREAD_COUNT in 10
do
    export OMP_NUM_THREADS=$THREAD_COUNT
    export FNAME= ${COUNTER}

    ../main -r $RADIUS -S $SPLIT_RATIO -o out.${counter}.json $FILENAME
    counter=$((counter+1))
done

python combine_outputs.py > results.json

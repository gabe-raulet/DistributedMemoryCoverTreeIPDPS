#!/bin/bash

# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads

export FILENAME=../points
export RADIUS=0.05
export SPLIT_RATIO=0.5

counter=1

for THREAD_COUNT in 1 2 4 8
do
    export OMP_NUM_THREADS=$THREAD_COUNT
    export FNAME= ${COUNTER}

    ../main -r $RADIUS -S $SPLIT_RATIO -o out.${counter}.json $FILENAME
    counter=$((counter+1))
done

python combine_outputs.py > results.json

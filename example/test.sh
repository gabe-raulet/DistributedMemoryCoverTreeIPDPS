#!/bin/bash

# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads

export FILENAME=../../points
export RADIUS=2.0
export SPLIT_RATIO=0.8

counter=1

for THREAD_COUNT in 1 2 3 4 5 6 7 8 9 10 11 12
do
    export OMP_NUM_THREADS=$THREAD_COUNT
    export FNAME= ${COUNTER}

    ../../main -r $RADIUS -S $SPLIT_RATIO -o out.${counter}.json $FILENAME
    counter=$((counter+1))
done

python combine_outputs.py > results.json

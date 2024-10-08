#!/bin/bash

# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads

# export FILENAME=sift_base_1000000.fvecs
# export RADIUS=0

# counter=1

# for SPLIT_RATIO in 0.65 0.7 0.75 0.8 0.85 0.9 0.95
# do
    # for THREAD_COUNT in 10
    # do
        # export OMP_NUM_THREADS=$THREAD_COUNT
        # export FNAME= ${COUNTER}

        # ../main -r $RADIUS -S $SPLIT_RATIO -o out.${counter}.json $FILENAME
        # counter=$((counter+1))
    # done
# done

# python combine_outputs.py > results.json

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

export RADIUS=1.15
export FILENAME=../points

counter=1

for SPLIT_RATIO in 0.7
do
    for THREAD_COUNT in 1 2 4 8 16 32 64 128 256
    do
        export OMP_NUM_THREADS=$THREAD_COUNT
        export FNAME= ${COUNTER}

        numactl -i0-7 ../main -r $RADIUS -S $SPLIT_RATIO -o out.${counter}.json $FILENAME
        counter=$((counter+1))
    done
done

python combine_outputs.py > results.json

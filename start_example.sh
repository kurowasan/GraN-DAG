#!/bin/bash
# make sure container.simg is in your working directory
# make sure singularity is installed on your computer.

CODE_PATH="/your/path/to/code/"  # path to the code root
EXP_PATH="/your/path/to/experiments/"  # this folder will contain all artifacts saved by the program
DATA_PATH="/your/path/to/data/"  # this should contain (data1.npy, DAG1.npy, CPDAG1.npy), (data2.npy, DAG2.npy, CPDAG2.npy), ...
MODEL="NonLinGaussANM"  # or NonLinGauss
DATA_INDEX=1  # Choose which dataset to use. Program will use (data${DATA_INDEX}.npy, DAG${DATA_INDEX}.npy, CPDAG${DATA_INDEX}.npy)
NUM_VARS=10  # should match the data provided

# Run python program. This is the Augmented Lagrangian procedure, without PNS. (add --pns to the command line to perform PNS)
singularity exec --containall -B $DATA_PATH:/dataset/ -B $EXP_PATH:/final_log/ -B $CODE_PATH:/code/ ./container.simg bash -c "cd /code && python main.py --exp-path /final_log/ --data-path /dataset/ --i-dataset ${DATA_INDEX} --model $MODEL --train --to-dag --num-vars ${NUM_VARS} --jac-thresh"

# To use the baselines methods, replace main.py by the right path (e.g. for notears, it is notears/main.py)
# singularity exec --containall -B $DATA_PATH:/dataset/ -B $EXP_PATH:/final_log/ -B $CODE_PATH:/code/ ./container.simg bash -c "cd /code && python cam/main.py --exp-path /final_log/ --data-path /dataset/ --i-dataset ${DATA_INDEX}"

#!/bin/bash
# TODO: finalize...
# make sure <TODO: insert container name> is in your working directory
# make sure singularity is installed on your computer.

CODE_PATH="/your/path/to/code/"  # path to the code root
EXP_PATH="/your/path/to/experiments/"  # this folder will contain all artifacts saved by the program
DATA_PATH="/your/path/to/data/"  # this should contain (data1.npy, DAG1.npy, CPDAG1.npy), (data2.npy, DAG2.npy, CPDAG2.npy), ...
MODEL="NonLinGaussANM"  # or NonLinGauss
DATA_INDEX=1  # Choose which dataset to use. Program will use (data${DATA_INDEX}.npy, DAG${DATA_INDEX}.npy, CPDAG${DATA_INDEX}.npy)
NUM_VARS=10  # should match the data provided

# Run python program. This is the Augmented Lagrangian procedure, without PNS. (add --pns to the command line to perform PNS)
# TODO: watch container name + paths inside the container probably don't exist.
singularity exec --containall -B $DATA_PATH:/dataset/ -B $EXP_PATH:/final_log/ -B $CODE_PATH:/code/ ./pytorch_r.sif bash -c "/opt/miniconda3/envs/pytorch/bin/python /code/main.py --exp-path /final_log/ --data-path /dataset/ --i-dataset ${DATA_INDEX} --model $MODEL --train --to-dag --num-vars ${NUM_VARS} --jac-thresh"
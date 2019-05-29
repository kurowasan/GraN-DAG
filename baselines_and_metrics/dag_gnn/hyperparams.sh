#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

module load python/3.6
BASEDIR="$PWD"

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index numpy scipy tqdm networkx torch

cd $BASEDIR
mkdir $SLURM_TMPDIR/data $SLURM_TMPDIR/results

DATASET="$(basename "$1" .zip)"
DAG="hyperparams_${DATASET}_${SLURM_ARRAY_TASK_ID}"

cp $1 $SLURM_TMPDIR/data
unzip $SLURM_TMPDIR/data/$DATASET.zip -d $SLURM_TMPDIR/data
mkdir $SLURM_TMPDIR/results/$DAG

python train.py $SLURM_TMPDIR/data/$DATASET --file_id 0 --save-folder $SLURM_TMPDIR/results/$DAG --encoder-hidden $SLURM_ARRAY_TASK_ID --decoder-hidden $SLURM_ARRAY_TASK_ID

cp -r $SLURM_TMPDIR/results/$DAG $SCRATCH/causal-learning/results_dag_gnn

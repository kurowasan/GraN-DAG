# Gradient-Based Neural DAG Learning

This code was written by the authors of the NeurIPS 2019 submission: "Gradient-Based Neural DAG Learning". Our implementation is in PyToroch. The directory `/baselines_and_metrics` contains the code we used for our baselines. DAG-GNN and NOTEARS are both in python (DAG-GNN uses PyTorch). All other baselines are in R programming language. Note that the metrics we report are also implemented in R. All R code is in `/baselines_and_metrics/rcode`.

## Run the code
To test GraN-DAG, you will only need to install Singularity to use the container image provided (pytorch_r.sif). TODO: provide link to installation instruction. Once this is done, you can start experimenting easily with the file `start_example.sh`.

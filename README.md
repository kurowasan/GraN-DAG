# Gradient-Based Neural DAG Learning

This code was written by the authors of the NeurIPS 2019 submission: "Gradient-Based Neural DAG Learning". Our implementation is in PyTorch. The directory `/baselines_and_metrics` contains the code we used for our baselines. DAG-GNN and NOTEARS are both in python (DAG-GNN uses PyTorch). All other baselines are in R programming language. Note that the metrics we report are also implemented in R. All R code is in `/baselines_and_metrics/rcode`.

## Run the code
To test GraN-DAG, you only need to install Singularity (version >= 3.0) to use the container image pytorch_r.sif. To install singularity, use the following instructions: [https://www.sylabs.io/guides/3.0/user-guide/installation.html](https://www.sylabs.io/guides/3.0/user-guide/installation.html). To download our container image pytorch_r.sif, use this link: [https://gofile.io/?c=BJvCP2](https://gofile.io/?c=BJvCP2). Once this is done, you can start experimenting easily with the file `start_example.sh` (make sure to change the paths).

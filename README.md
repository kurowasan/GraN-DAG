# Gradient-Based Neural DAG Learning

This code was written by the authors of the ICLR 2020 submission: "Gradient-Based Neural DAG Learning". 
Our implementation is in PyTorch but some functions rely on the 
[Causal Discovery Toolbox](https://diviyan-kalainathan.github.io/CausalDiscoveryToolbox/html/index.html) which relies 
partly on the R programming language. 

## Run the code
To use our implementation of GraN-DAG, simply install Singularity (instructions: [https://www.sylabs.io/guides/3.0/user-guide/installation.html](https://www.sylabs.io/guides/3.0/user-guide/installation.html)) 
and run the code in our container (download it here: [https://zenodo.org/record/7700409](https://zenodo.org/record/7700409)). Use start_example.sh (update the paths) to launch the differents methods (GraN-DAG, DAG-GNN, NOTEARS, CAM).

MPI-SGD
=======

Example usage:

1. Put your data into a directory called 'data' in the root of the repository (or change the paths in the file 'scripts/run_single_experiment.sh' accordingly
2. Run 'make' in the root directory
3. Run 'bash scripts/run_single_experiment.sh higgs all_reduce decreasing sparse 4 1 20 1' in the root directory (ommit the arguments to see all possibilities)

- For running multiple experiments and store the results please make use of the bash script 'scripts/run_experiments_skeleton.sh' or have a look at 'scripts/run_all_combinations.sh
- The python script 'scripts/produce_graphs.py' provides some functions to generate graphics based on the multiple experiments run using this skeleton script

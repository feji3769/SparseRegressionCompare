# SparseRegressionCompare
Comparing sparse bayesian regression methods. 

# Dependencies : 

1. [argparse](https://cran.r-project.org/web/packages/argparse/readme/README.html)

2. [varbvs](https://cran.r-project.org/web/packages/varbvs/index.html)


The horseshoe code is simply a copy of the code from [horseshoe](https://cran.r-project.org/web/packages/horseshoe/index.html) on Cran. 

# Run experiments : 
## Single Run
Rscript Experiments.R -n 100 -p 200 -SNR 4 -method MCMCHorseshoe


This will create and save results to "experimentResults/experiment_results.csv", or append if it already exists. 

## Batch of Experiments
running the bash script runExperiment will run 50 replicates of 8 different settings (changing n, p, and SNR) for one model. Changing the model argument in this bash script will run the same experiment for that inference technique. 

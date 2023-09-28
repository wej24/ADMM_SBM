# Simultaneous estimation of connectivity and dimensionality in samples of networks
This repository contains R code to run the ADMM algorithm and simulations in the "Simultaneous estimation of connectivity and dimensionality in Samples of Networks" paper. 

ADMM.R file in the Algorithm folder contains the main ADMM algorithm with a single $\lambda$ value.

cvADMM.R file in the Algorithm folder contains the M-fold cross-validation of the ADMM algorithm for multiple adjacency matrices with a sequence of $\lambda$ values. 

simulationFunc.R and simulationFuncHeSBM contain the helper functions to run the simulations in the paper. 

R files in the HoSBM folder and HeSBM folder contain the code to reproduce the simulations in the paper. To reproduce the simulations, simply run the R file in these two folders. 

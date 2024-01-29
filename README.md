# Simultaneous estimation of connectivity and dimensionality in samples of networks
This repository contains R code to run the ADMM algorithm and simulations in the "Simultaneous estimation of connectivity and dimensionality in Samples of Networks" paper. 

ADMM.R file in the Algorithm folder contains the main ADMM algorithm with a single $\lambda$ value.

cvADMM.R file in the Algorithm folder contains the M-fold cross-validation of the ADMM algorithm for multiple adjacency matrices with a sequence of $\lambda$ values. 

simulationFunc.R and simulationFuncHeSBM contain the helper functions to run the simulations in the paper. 

R files in the HoSBM folder and HeSBM folder contain the code to reproduce the simulations in the paper. To reproduce the simulations, simply run the R file in these two folders. 

real_data.R file in the real_data folder contains the real data analysis. If researchers are interested in running this file with the data used in the paper, you can request the data by contacting the authors of "Bias-adjusted spectral clustering in multi-layer stochastic block models. Journal of the American Statistical Association 0(0), 1–13" by Jing Lei and Kevin Z. Lin. The Bias-adjusted spectral clustering code can be found at https://github.com/linnykos/networkSoSD. 

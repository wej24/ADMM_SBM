# Simultaneous estimation of connectivity and dimensionality in samples of networks
This repository contains R code to run the ADMM algorithm and simulations in the paper "Simultaneous estimation of connectivity and dimensionality in Samples of Networks" by Wenlong Jiang, Chris McKennan, Jesús Arroyo, and Joshua Cape

# Example for how to use algorithm
algorithm_example.R shows a simple example for how to use our algorithm

# Simulations in the paper
Note: The code only run one replicate and one case for each setting in the paper since we run all results on the cluster and each replicate and each case run independently. 

1. Figure1.R contains the code to plot Figure 1 in the main text.
2. Figure7.R contains the code to plot Figure 7 in the supplement.
3. MonoSBM_K10d1_TrueZ.R contains the code to run the simulation in the section 5.1 in the main text.
4. MonoSBM_K10d2_estZ.R contains the code to run the simulation in the section 5.2 in the main text.
5. MonoSBM_K10d2_reEstZ.R contains the code to run the simulation in the section 5.3 in the main text.
6. MultiSBM_K3_estZ.R contains the code to run the simulation in the section 5.4 in the main text.
7. real_data.R contains the code to run the real data analysis in the section 6 in the main text.
8. MonoSBM_single_K2d1_TrueZ.R and MonoSBM_single_K2d2_TrueZ.R contain the code to run the simulation in the section S3.1 in the supplement.
9. MonoSBM_K2d1_TrueZ.R and MonoSBM_K2d2_TrueZ.R contain the code to run the simulation in the section S3.2 in the supplement.
10. MonoSBM_K2d1_estZ.R and MonoSBM_K2d2_estZ.R contain the code to run the simulation in the section S3.3 in the supplement.

# More on real data
If researchers are interested in the original real data used in the paper, you can request the data by contacting the authors of "Bias-adjusted spectral clustering in multi-layer stochastic block models. Journal of the American Statistical Association 0(0), 1–13" by Jing Lei and Kevin Z. Lin. The Bias-adjusted spectral clustering code can be found at https://github.com/linnykos/networkSoSD. 

# Simultaneous estimation of connectivity and dimensionality in samples of networks
This repository contains R code for the paper "Simultaneous estimation of connectivity and dimensionality in samples of networks" by Wenlong Jiang, Chris McKennan, Jesús Arroyo, and Joshua Cape.

# Example for how to use the ADMM algorithm
algorithm_example.R provides a simple example for how to use the algorithm.

# Simulations in the paper
Note: The code only runs one replicate and one case for each setting in the paper. For the paper, we used the University of Pittsburgh computing cluster, with each replicate and each case run independently. 

1. Figure1.R contains the code to plot Figure 1 in the main text.
2. Figure7.R contains the code to plot Figure 7 in the supplement.
3. MonoSBM_K10d1_TrueZ.R contains the code to run the simulation in Section 5.1 in the main text.
4. MonoSBM_K10d2_estZ.R contains the code to run the simulation in Section 5.2 in the main text.
5. MonoSBM_K10d2_reEstZ.R contains the code to run the simulation in Section 5.3 in the main text.
6. MultiSBM_K3_estZ.R contains the code to run the simulation in Section 5.4 in the main text.
7. MonoSBM_single_K2d1_TrueZ.R and MonoSBM_single_K2d2_TrueZ.R contain the code to run the simulation in Section S3.1 in the supplement.
8. MonoSBM_K2d1_TrueZ.R and MonoSBM_K2d2_TrueZ.R contain the code to run the simulation in Section S3.2 in the supplement.
9. MonoSBM_K2d1_estZ.R and MonoSBM_K2d2_estZ.R contain the code to run the simulation in Section S3.3 in the supplement.

# Real data in the paper
real_data.R contains the code to run the real data analysis in Section 6 in the main text.

The original primate brain dataset can be obtained by contacting the authors of "Bias-adjusted spectral clustering in multi-layer stochastic block models" (JASA 2023), Jing Lei and Kevin Z. Lin. The corresponding bias-adjusted spectral clustering code can be found at https://github.com/linnykos/networkSoSD.

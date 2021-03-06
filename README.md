# sdbmpABC
Sample code for the paper by Buckwar, Tamborrino &amp; Tubikanec "Spectral Density-Based and Measure-Preserving ABC for partially observed diffusion processes An illustration on Hamiltonian SDEs", Stat. Comput. 30, 627-648, 2020, https://link.springer.com/content/pdf/10.1007/s11222-019-09909-6.pdf

The code was written by Irene Tubikanec (firstname dot secondname at jku.at) and then updated by Massimiliano Tamborrino (firstname dot secondname at jku.at).

# What can you find in the package
Since exact simulation for the models of Section 4 is available, here we provide the code for performing the Spectral Density-Based and Measure-Preserving ABC (sdbmpABC) for the Jansen and Rit Neural Mass Model (JR-NMM) (25) of Section 5. In particular, we provide what we denoted Algorithm 1 (ii) in the paper to reproduce Figure 8. 

The proposed acceptance-rejection ABC algorithm is based on two key ingredients: 
1) summary statistics based on the estimated invariant density and invariant spectral density; 
2) a structure-preserving numerical Strang splitting method. 
The numerical method preserves the measure properties of the model, guaranteeing the succesfull inference.

# How to install the package
The simplest way is to do it via devtools, using devtools::install_github("massimilianotamborrino/sdbmpABC")

# How does the package work
These are the steps needed to reproduce the inference of θ = (σ,μ,C) of the JR-NMM based on simulated reference data, i.e., to  reproduce Figure 8 (with pairwise scatterplots of the kept ABC posterior samples in addition).

1. Install the provided package "sdbmpABC" (see above).
2. Specify the ABC setting (cut, N, M, w, T, h) in the file "run_ABC_JRNMM_sigmuC.R".
   The pre-defined setting  
   cut=10^3 (corresponding to the kept posterior samples)
   N=2.5*10^6 (corresponding to the number of total samples from the uniform priors)
   M=30 (corresponding to the number of observed datasets, i.e., number of trajectories of the output process.
   w=1930.17 (weight entering into the distance)
   T=200 (length of the interval in which to simulate the data)  
   h=0.002 (time step)
   is the one used in the corresponding manuscript.
4. Run the file "run_ABC_JRNMM_sigmuC.R".
5. Check if the kept (marginal) posterior samples are successfully stored 
   as txt files (sig, mu, C) in the working folder.
6. Run the file "visualisation_of_results.R".
7. The figure visalising the ABC marginal posterior densities of θ = (σ,μ,C) (Top panels) and the pairwise scatterplots of the kept posterior samples (Lower panels) of the stochastic JR-NMM (25) obtained from Algorithm 1 (ii) is now stored. The horizontal red lines and the vertical black lines represent the uniform priors and the true parameter values, respectively. Note that the Top panels correspond to Figure (8) of the manuscript.

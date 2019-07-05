# sdbmpABC
Sample code for the paper by Buckwar, Tamborrino &amp; Tubikanec (2019) "Spectral Density-Based and Measure-Preserving ABC for partially observed diffusion processes An illustration on Hamiltonian SDEs", Preprint at https://arxiv.org/abs/1903.01138

The code was written by Irene Tubikanec (firstname dot secondname at jku.at) and then updated by Massimiliano Tamborrino (firstname dot secondname at jku.at).

Since exact simulation for the models of Section 4 is available, here we provide the code for performing the Spectral Density-Based and Measure-Preserving ABC (sdbmsABC) for the Jansen and Rit Neural Mass Model (JR-NMM) (25) of Section 5.

These are the steps needed to reproduce the inference of theta=(sigma,mu,C) of the JR-NMM based on simulated reference data:

1. Install the provided package "Splitting".
2. Install the packages called in the file "ABC.R".
3. Specify the ABC setting (cut, N, M, w, T, h) in the file "run_ABC_JRNMM_sigmuC.R".
   The pre-defined setting 
   cut=10^3
   N=2.5*10^6
   M=30
   w=1930.17
   T=200
   h=0.002
   is the one used in the corresponding manuscript.
4. Run the file "run_ABC_JRNMM_sigmuC.R".
5. Check if the kept (marginal) posterior samples are successfully stored 
   as txt files (sig, mu, C) in the folder "Marginal_posteriors".
6. Run the file "visualisation_of_results.R".
7. The figure visalising the results is now stored. 


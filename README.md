# Markov-Jump-Adjusted-Particle-Filter-MJAPF

This implementation of a novel Particle Filter is based on the work presented in the following paper:

Hedayati M., Rahimi A. (2025) A Hybrid Framework for Real-Time Satellite Fault Diagnosis Using Markov Jump-Adjusted Models and 1D Sliding Window Residual Networks, Acta Astronautica, from Elsevier BV, Volume 228, March 2025, Issue N/A, pp. 1066-1087, doi: 10.1016/j.actaastro.2024.12.057

For a better understanding of how this state estimation algorithm works, please consult the paper.

In this repository, a double-state estimator implementation of the Markov Jump-Adjusted Particle Filter (MJAPF) algorithm is provided.

The 'pf_m.m' script contains the function for calling the MJAPF algorithm. 'model_output.m' is the script responsible for propagating a dynamic system's model based on an individual particle trajectory up to that time-step. 'jumping_strategy.m' contains the code on deciding the jumps based on the transition probability matrix. The 'process_model_kde.m' script obtains the approximated Process noise PDF of a dynamic system based on Kernel Density Estimation (KDE) given that a large number of dynamic system simulations are possible. 

Parallel Pooling is needed to run this code on MATLAB for efficiency, therefore, the Parallel Computing License is needed for that. If you do not have it, you can revert the "parfor"s to standard for-loops.


If you have any issues with the codes or questions, please contact me at mohammadshedayati@gmail.com or hedayat2@uwindsor.ca

This project reuses the "subplotXmanyY_er.m" code from Dr. Afshin Rahimi's repository, available (https://github.com/drarahimi/MATLAB_Professional_Plots?tab=MIT-1-ov-file), licensed under the MIT License.

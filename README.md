# SARS-CoV-2-Evolution
Quantitative modelling of SARS-CoV-2 evolutionary dynamics.

## main.m
main.m runs one simulation with default input parameter values and plots the number of uninfected cells, infected cells and free viral particles over time.
A sparse sensitivity analysis is also run and plotted in main.m. Less parameter sets and replicates are included as running many simulations is probably not feasible on a single machine. Therefore, results will differ significantly compared to the PRCC matrix reported in the manuscript. 

## Sensitivity analysis
The Latin Hypercube Sampling initialization is done by running "Simulation_code/lhs_initialization.m". This script generates parameters sets which can be used as input parameter values to "Simulation_code/gillespie.m". The parameters which are varied include the infection rate *a*, the clearance/death rate *b*, the replication rate *r0* and the mutation rate *mu*. All fours parameters are varied 20% around their default values. Partial Rank Correlation Coefficient (PRCC) analysis can then be done within MATLAB with the *partialcorri* function.

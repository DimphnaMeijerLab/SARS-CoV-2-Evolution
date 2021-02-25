function run_gillespie_transmission_parallel(x,~)

    % init
    U0 = 1e4;
    mu_array = 1e-6;
    distribution = 'normal';
    sigma = 0.1;
    wholeGenome = false;
    fitnessFunction = @multiplicative_fitness;
    N_ind = 3;
    
    gillespie_transmission(str2double(x), U0, mu_array, N_ind, ...
                                    'distribution', distribution, ...
                                    'sigma', sigma, ...
                                    'fitnessFunction', fitnessFunction, ...
                                    'wholeGenome', wholeGenome);

end
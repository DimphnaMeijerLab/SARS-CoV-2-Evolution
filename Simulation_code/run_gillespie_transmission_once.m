function run_gillespie_transmission_once(nr)

    % init
    U0 = 1e5;
    mu = 1e-6;
    distribution = 'normal';
    sigma = 0.1;
    wholeGenome = false;
    fitnessFunction = @multiplicative_fitness;
    downloadDate = '20210120';
    mistmatchThreshold = 't3';
    
    gillespie_transmission(str2double(nr), U0, mu, ...
                           downloadDate, mistmatchThreshold,...
                           'distribution', distribution,...
                           'sigma', sigma, ...
                           'fitnessFunction', fitnessFunction, ...
                           'wholeGenome', wholeGenome);
                                        
end
function run_gillespie_once(nr)
    
    % init
    U0 = 1e4;
    mu_array = logspace(-8,-3,20);
    distribution = 'normal';
    wholeGenome = false;
    sigma = 0.1;
    fitnessFunction = @multiplicative_fitness;
    downloadDate = '20210120';
    mistmatchThreshold = 't3';
    
    [data, params] = gillespie(str2double(nr), U0, mu_array, ...
                                downloadDate, mistmatchThreshold,...
                                'distribution', distribution, ...
                                'sigma', sigma, ...
                                'fitnessFunction', fitnessFunction, ...
                                'wholeGenome', wholeGenome);
                                
    
 	fname = ['Output/Gillespie', nr, '.mat'];
    save(fname, 'data', 'params', '-v7.3');
                                        
end
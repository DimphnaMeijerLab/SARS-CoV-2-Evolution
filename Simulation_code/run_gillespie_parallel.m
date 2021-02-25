function run_gillespie_parallel(x1, x2_max) 

    x2_max = str2double(x2_max);

    % init
    U0 = 1e4;
    mu_array = 1e-6;
    distribution = 'normal';
    sigma = 0.1;
    wholeGenome = false;
    fitnessFunction = @multiplicative_fitness;
    
    for x2 = 0:x2_max
         nr = [x1, num2str(x2)];
         fprintf(['\n Starting iteration nr -> ',nr,' \n'])
         % run gillespie
         [data, params] = gillespie(str2double(nr), U0, mu_array, ...
                                    'distribution', distribution, ...
                                    'sigma', sigma, ...
                                    'fitnessFunction', fitnessFunction, ...
                                    'wholeGenome', wholeGenome);
         % save
         fname = ['Output/Gillespie', nr, '.mat'];
         save(fname, 'data', 'params', '-v7.3');
    end
end
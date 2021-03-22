function run_gillespie_parallel(x1, x2_max) 

    % mkdir('./Output');
    x2_max = str2double(x2_max);
    
    % init
    U0 = 1e6;
    mu_array = 1e-6;
    distribution = 'normal';
    sigma = 0.1;
    wholeGenome = false;
    fitnessFunction = @multiplicative_fitness;
    downloadDate = '20210120';
    mistmatchThreshold = 't3';
    
    % Shut down par-pool if it still existed
    % delete(gcp('nocreate'))
    % Open a new parpool
    % parpool(8)
    
    tStart = tic;
    for x2 = 0:x2_max
        
        nr = [x1, num2str(x2)];
        fprintf(['\n Starting iteration nr -> ',nr,' \n'])
        
        % run gillespie
        [data, params] = gillespie(str2double(nr), U0, mu_array, ...
                                   downloadDate, mistmatchThreshold,...
                                   'distribution', distribution, ...
                                   'sigma', sigma, ...
                                   'fitnessFunction', fitnessFunction, ...
                                   'wholeGenome', wholeGenome);
                                
        fname = ['Output/Gillespie', nr, '.mat'];
        parsave(fname, data, params);
    end
    
    fprintf(['\n Total simulation time: ',num2str(toc(tStart)),' \n'])
    % Shut down par-pool
    delete(gcp('nocreate'))
    
end

function parsave(fname, data, params)
    save(fname, 'data', 'params', '-v7.3');
end
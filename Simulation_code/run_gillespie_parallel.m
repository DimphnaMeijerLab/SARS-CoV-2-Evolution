function run_gillespie_parallel(x1) 
    % init
    U0 = 1e4;
    mu_array = 1e-6;
    
    for x2 = 0:4
         nr = [x1, num2str(x2)];
         % run gillespie
         [data, params] = gillespie(str2double(nr), U0, mu_array);
         % save
         fname = ['Output/Gillespie', nr, '.mat'];
         save(fname, 'data', 'params', '-v7.3');
    end
end
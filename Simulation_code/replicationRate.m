function r = replicationRate(d, r0, distribution, fitnessFunction, sigma, beta)
    %----------------------------------------------------------------------
    %
    % This function determines the fitness of a strain with a distance d from
    % the reference sequence, which has fitness r0.
    % The fitness can be chosen to be (1) a normal distribution with 
    % mean mu (predicted by logistic regression) and standard deviation
    % sigma, or (2) a gamma distribution with mean mu and standard
    % deviation sigma.
    %
    % INPUT PARAMETERS
    % ----------------
    %
    % d: nProtx1 array
    % Hamming distance of strain i for all proteins.
    %
    % r0: double
    % Replication rate of reference sequence.
    %
    % distribution: string; choose from 'normal' or 'gamma' or 'empirical'.
    % Probability distribution of the fitness of newly emerged strains.
    %
    % fitnessFunction: function handle
    % Function to calculate fitness. The first input argument must be b (a
    % 2x1 column vector with b0 and b1, the outputs of the logistc fit).
    % The second input argument must be the hamming distance d.
    %
    % sigma: array or double
    % Standard deviation of normal distribution (if distribution is
    % 'normal') or standard deviation of gamma distribution (if
    % distribution is 'gamma'). Enter an array to specify sigma per protein
    % or a single float to use the same sigma for each protein.
    %    
    % beta (2xnProt array, where nProt is the number of proteins).
    % The first row contains beta0 for all proteins, the second row
    % contains beta1 for all proteins.
    %
    % OUTPUT PARAMETER
    % ----------------
    %
    % r: double
    % Replication rate of the strain.
    %
    %----------------------------------------------------------------------
    N = length(d);
    if length(sigma) == 1
        sigma = zeros(1,N) + sigma;
    end
    
    rArray = zeros(N,1);
    for i=1:N
        if d(i) == 0
            replRate = 1;
        else
            b = beta(:, i);
            mu = fitnessFunction(b, d(i));
            
            if strcmp(distribution, 'normal') 
                pd = makedist('Normal','mu',mu,'sigma',sigma(i));
                t = truncate(pd,0,inf);
                replRate = random(t); 
                
            elseif strcmp(distribution, 'gamma')     
                shape = mu^2 / sigma(i);
                scale = sigma(i) / mu;
                replRate = random('Gamma', shape, scale);
            
            elseif strcmp(distribution, 'empirical')
                % Make probability distributions and truncate them
                if rand(1) < 0.958      % the probability of a deleterious mutation is 95.8%
                    mean = mu * 0.1904; % 0.1904 is the expected fitness reduction of Sanjuan. We downscale it according to the expected fitness
                    mu_param = log( mean ) - 1.206^2 / 2;
                    pd_deleterious = makedist('Lognormal','mu',mu_param,'sigma',1.206);
                    pd_deleterious = truncate(pd_deleterious,0,mu);
                    replRate = random(pd_deleterious);
                else
                    pd_beneficial = makedist('Exponential','mu',mu);
                    pd_beneficial = truncate(pd_beneficial,mu,inf);
                    replRate = random(pd_beneficial);
                end

            else
                error('Invalid distribution name');
            end
        end
        
        rArray(i) = replRate;
    end
    r = r0 * prod(rArray);
    %r = rArray;
end
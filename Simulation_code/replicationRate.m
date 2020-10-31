function r = replicationRate(d, r0, sigma, beta, distribution)
    %----------------------------------------------------------------------
    
    % This function determines the fitness of a strain with a distance d from
    % the reference sequence, which has fitness r0.
    % The fitness can be chosen to be (1) a normal distribution with 
    % mean mu (predicted by logistic regression) and standard deviation
    % sigma, or (2) a gamma distribution with mean mu and standard
    % deviation sigma.
    
    % INPUT PARAMETERS
    % ----------------
    
    % d: nProtx1 array
    % Hamming distance of strain i for all proteins.
    
    % r0: double
    % Replication rate of reference sequence.
    
    % beta: 2xnProt array, where nProt is the number of proteins.
    % The first row contains beta0 for all proteins, the second row
    % contains beta1 for all proteins.
    
    % distribution: string; choose from 'normal' or 'gamma'.
    % Probability distribution of the fitness of newly emerged strains.
    
    % OUTPUT PARAMETER
    % ----------------
    
    % r: double
    % Replication rate of the strain.
    
    %----------------------------------------------------------------------
    N = length(sigma);
    rArray = zeros(N,1);
    for i=1:N
        if d(i) == 0
            replRate = 1;
        else
            b = beta(:, i);
            mu = predict(b, d(i));
            
            if strcmp(distribution, 'normal') 
                pd = makedist('Normal','mu',mu,'sigma',sigma(i));
                t = truncate(pd,0,inf);
                replRate = random(t); 
                
            elseif strcmp(distribution, 'gamma')     
                shape = mu^2 / sigma(i);
                scale = sigma(i) / mu;
                replRate = random('Gamma', shape, scale);
                
            else
                disp('Invalid distribution name');
            end
        end
        
        rArray(i) = replRate;
    end
    r = r0 * prod(rArray);
    %r = rArray;
end
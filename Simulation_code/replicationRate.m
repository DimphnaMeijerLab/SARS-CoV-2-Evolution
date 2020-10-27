function r = replicationRate(d, r0, sigma, beta, distribution)
% This function determines the fitness of a strain with a distance d from
% the reference sequence, which has fitness r0.
% The fitness is normal distributed with mean mu (predicted by logistic 
% regression) and standard deviation sigma.
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
function w = multiplicative_fitness(b, d)
    d0 = - b(1) ./ b(2);            % d0 = -b0/b1
    c = 2 ^ (-1/d0);                % choose c s.t. fitness(d0) = 1/2    
    w = c .^ d;
end

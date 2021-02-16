function w = truncation_model_fitness(b, d)
    d0 = - b(1) ./ b(2);              % d0 = -b0/b1
    w = zeros(size(d));
    w(d <= d0) = 1;
    w(d > d0) = 0;
end
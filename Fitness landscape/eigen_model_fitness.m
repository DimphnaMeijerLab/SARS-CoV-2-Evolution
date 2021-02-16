function w = eigen_model_fitness(b, d)
    w = zeros(size(d));
    w(d==0) = 1;
    w(d>0) = 0.5;
end
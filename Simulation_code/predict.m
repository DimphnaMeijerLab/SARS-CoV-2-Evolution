function p = predict(b, X)
    m = size(X, 1);
    X = [ones(m, 1) X];
    
    p = sigmoid(X * b);
end

function g = sigmoid(z)
    g = 1 ./ (1 + exp(-z));
end

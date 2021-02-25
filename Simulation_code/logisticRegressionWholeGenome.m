function [b,sigma] = logisticRegressionWholeGenome(mismatchBooleanWholeGenome, lambda)
    y = mismatchBooleanWholeGenome;
    N = length(y);
    X = transpose(0 : N-1);
    % Do logistic regression
    Mdl = fitclinear(X,y,'Learner','logistic','BetaTolerance',1e-10,'Lambda',lambda);
    b = [Mdl.Bias; Mdl.Beta];
    prediction = predict(b, X);
   	sigma = sqrt (sum((prediction - y) .^2) / (N - 2) );
end
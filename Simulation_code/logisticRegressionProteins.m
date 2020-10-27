function [beta, sigma] = logisticRegressionProteins()
% Function performs a logistic regression on the number of mutations for 
% protein. The mutation numbers are stored in 'mismatchBoolean'.
% beta0 and beta1 are arrays with the values of beta1 and beta2 for all
% proteins.
% Sigma is an array with the standard error of the logistic regression for
% each protein

    % load the data
    data = load('mismatchBoolean0923.mat');
    mismatchBoolean = data.mismatchBoolean;
    proteinNames = fields(mismatchBoolean);
    lambda = 0.001;
    % initialize beta arrays:
    beta = zeros(2,length(proteinNames));
    sigma = zeros(1, length(proteinNames));
    % loop over proteins
    for i = 1:length(proteinNames)
        protein = proteinNames{i};
        y = mismatchBoolean.(protein)';
        N = length(y);
        X = transpose(0 : N-1);
 
        % Do logistic regression
        %b = logisticRegression(X, y);
        %b = lassoglm(X,y,'binomial','link','logit','Lambda', logspace(-15,0,100));
        Mdl = fitclinear(X,y,'Learner','logistic','Regularization','ridge','Lambda',lambda);
        b = [Mdl.Bias; Mdl.Beta];
        beta(:,i) = b;
        % Find standard error
        prediction = predict(b, X);
        sigma(i) = sqrt (sum((prediction - y) .^2) / (N - 2) );
    end
end


clear all
close all

downloadDate = '20210120';

addpath('../Data')
fileName = ['mismatchBooleanOverTime',downloadDate,'_t3.mat'];
load(fileName)

mismatchBoolean = mismatchBooleanOverTime(end);
proteinNames = fields(mismatchBoolean);
numProteins = length(proteinNames);

lambda = 0;
[betaFit, sigmaFit] = logisticRegressionProteins(mismatchBoolean, lambda);

d0Fit = - betaFit(1,:) ./ betaFit(2,:);    % d0 = -b0/b1
alfaFit = - betaFit(2,:);                     % alfa = -b1

y = mismatchBoolean.NSP1';
N = length(y);
X = transpose(0 : N-1);
s1 = calculate_sigma(alfaFit(1), d0Fit(1), X, y);
s2 = calculate_sigma(0, d0Fit(1), X, y);

function [beta, sigma] = logisticRegressionProteins(mismatchBoolean, lambda)
    
    %----------------------------------------------------------------------
    % Function performs a logistic regression on the number of mutations for 
    % each protein. The mutation numbers are stored in 'mismatchBoolean'.
    %
    % The assumed relationship between expected fitness contribution <f> 
    % and hamming distance d is
    % <f> = 1 ./ (1 + exp(-(beta0 + beta1*d)));
    % This function outputs the fitted beta0 and beta1 for each protein and
    % the least-square errors of the fit.
    %
    % INPUT PARAMETERS
    % ----------------
    %
    %   mismatchBoolean (struct)
    %   Data structure array with the mismatch booleans of all
    %   proteins.
    %
    %   lambda (double)
    %   L2-norm (ridge) regression parameter.
    %
    % OUTPUT PARAMETERS
    % ----------------
    %
    %   beta (2xnProt array, where nProt is the number of proteins).
    %   The first row contains beta0 for all proteins, the second row
    %   contains beta1 for all proteins.
    
    %   sigma (1xnProt array)
    %   An array with the standard error of the logistic regression for
    %   each protein.
    %----------------------------------------------------------------------

    % load the data
    proteinNames = fields(mismatchBoolean);
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
        Mdl = fitclinear(X,y,'Learner','logistic','Regularization','ridge','Lambda',lambda);
        b = [Mdl.Bias; Mdl.Beta];
        beta(:,i) = b;
        % Find standard error
        prediction = predict(b, X);
        sigma(i) = sqrt (sum((prediction - y) .^2) / (N - 2) );
    end
end

%--------------------------------------------------------------------------
% Helper functions for logistic regression.
%--------------------------------------------------------------------------

function p = predict(b, X)
    m = size(X, 1);
    X = [ones(m, 1) X];
    
    p = sigmoid(X * b);
end

function sigma = calculate_sigma(alfa, d0, X, y)
    N = length(y);
    p = sigmoid(alfa*(d0 - X));
    sigma = sqrt (sum((p - y) .^2) / (N - 2) );
end

function g = sigmoid(z)
    g = 1 ./ (1 + exp(-z));
end

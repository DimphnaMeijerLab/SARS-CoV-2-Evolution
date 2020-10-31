function [beta, sigma] = logisticRegressionProteins()

    %----------------------------------------------------------------------
    % Function performs a logistic regression on the number of mutations for 
    % each protein. The mutation numbers are stored in 'mismatchBoolean'.
    
    % The assumed relationship between expected fitness contribution <f> 
    % and hamming distance d is
    % <f> = 1 ./ (1 + exp(-(beta0 + beta1*d)));
    % This function outputs the fitted beta0 and beta1 for each protein and
    % the least-square errors of the fit.
    
    % OUTPUT PARAMETERS
    % ----------------
    
    % beta: 2xnProt array, where nProt is the number of proteins.
    % The first row contains beta0 for all proteins, the second row
    % contains beta1 for all proteins.
    
    % sigma: 1xnProt array
    % An array with the standard error of the logistic regression for
    % each protein.
    %----------------------------------------------------------------------

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
        Mdl = fitclinear(X,y,'Learner','logistic','Regularization','ridge','Lambda',lambda);
        b = [Mdl.Bias; Mdl.Beta];
        beta(:,i) = b;
        % Find standard error
        prediction = predict(b, X);
        sigma(i) = sqrt (sum((prediction - y) .^2) / (N - 2) );
    end
end


clear all
close all

downloadDate = '20210120';

addpath('../Data')
fileName = ['mismatchBooleanOverTime',downloadDate,'_t3.mat'];
load(fileName)

mismatchBoolean = mismatchBooleanOverTime(end);
proteinNames = fields(mismatchBoolean);
numProteins = length(proteinNames);

% Do logistic regression (without regularization)
lambda = 0;
[betaFit, sigmaFit] = logisticRegressionProteins(mismatchBoolean, lambda);
% Convert b0 and b1 into alfa and d0 parameters
d0Fit = - betaFit(1,:) ./ betaFit(2,:);         % d0 = -b0/b1
alfaFit = - betaFit(2,:);                       % alfa = -b1

s0 = 0.05;               % each protein should have sigma >= 0.1
tolerance = 1e-5;       % tolerance of root solving

for p = 1:length(proteinNames)
    protein = proteinNames{p};

    y = mismatchBoolean.(protein)';
    N = length(y);
    X = transpose(0 : N-1);

    % sigma = 0.1 should lie somewhere between alfa=0 and alfa=fitted_alfa.
    % If it doesn't, raise a message.
    d0 = d0Fit(p);
    s1 = calculate_sigma(alfaFit(p), d0, X, y);
    s2 = calculate_sigma(0, d0, X, y);
    
    fit_alfa = true;
    if (s1-s0) * (s2-s0) > 0
        disp(['Warning: Sigma is already larger than s0 in protein ', protein])
        fit_alfa = false;
    end
    % Find roots if sigma is smaller than 0.1
    if fit_alfa == true
        start_vector = [0,alfaFit(p)];  
        calculated_roots = find_root_with_bisection_method(start_vector,tolerance,s0,d0,X,y);
        alfaNew = calculated_roots(1);
    else % if not, use the alfa fitted by logistic regression
        alfaNew = alfaFit(p);
    end

    % Make plot
    maxy = find(y,1,'last') + 15;
    yax = y(1:maxy);
    Xax = linspace(0,length(yax)-1,length(yax))';
    xaxis = linspace(0,length(yax)-1, 10*length(yax))';
    b = [d0*alfaNew; -alfaNew];         % Convert alfa and d0 back in b0 and b1.
    S = predict(b, xaxis);
    S(1:2) = NaN;

    subplot(5,6,p)
    plot(xaxis,S, '-k')
    subplot(5,6,p)
    plot(Xax,yax,'.r','MarkerSize',10);
    hold on
    plot(xaxis, S, 'k','LineWidth',1);
    axis([0, max(Xax),0,1.1]);
    xticks([0,max(Xax)])   
    yticks([0,1])
    title(protein)
end

%%

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

%% Root solving with bisection method

function calculated_roots = find_root_with_bisection_method(start_vector,tolerance,s0,d0,X,y)
    %----------------------------------------------------------------------
    % This function finds the roots of the function sigma(alfa) - s0.
    % Sigma (error of logistic fit) is a function of alfa.
    % We want to find alfa s.t. sigma = s0. 
    %
    %   INPUTS
    %   ------
    %   start_vector (row vector with floats)
    %   Starting points for root finding. Supply 2 starting points for
    %   every root you want to find. Example: if you expect two roots, one
    %   around a=2 and one around a=4, then start_vector=[1.5 2.5 3.5 4.5].
    %
    %   tolerance (float)
    %   Root finding stops if the error is smaller than tolerance.
    %
    %   s0 (float)
    %   Desired error for each protein. We want to find alfa s.t. sigma = s0.
    %
    %   d0 (float)
    %   Original d0, found with logistic fit. This parameter does not
    %   change.
    %
    %   X (column vector with integers)
    %   Array with hamming distances.
    % 
    %   y (column vector, each value is either 0 or 1)
    %   Mismatch boolean of protein.
    %
    %   OUTPUTS
    %   -------
    %   calculated_roots (vector with floats)
    %   Found roots.
    %----------------------------------------------------------------------
    
    
    calculated_roots = zeros(1, length(start_vector)-1);

    for i = 1:length(start_vector)-1

        a1 = start_vector(i);
        a2 = start_vector(i+1);
        f1 = calculate_sigma(a1, d0, X, y) - s0;
        f2 = calculate_sigma(a2, d0, X, y) - s0;
        if f1 * f2 > 0
            exit('Root is not in between a1 and a2. Choose your inital a1 and a2 differently.')
        end

        error = 1;  % Initial error
        c = 0;      % counter

        while error > tolerance
            error = abs( (a2-a1)/(a1+a2) );
            a_sol = (a1 + a2) / 2;
            f_sol = calculate_sigma(a_sol, d0, X, y) - s0;
            f1 = calculate_sigma(a1, d0, X, y) - s0;
            f2 = calculate_sigma(a2, d0, X, y) - s0;

            % perform bisection
            if f1 * f_sol < 0      % solution is between a1 and a_sol
                a2 = a_sol;
            elseif f2 * f_sol < 0  % solution is between a2 and a_sol
                a1 = a_sol;
            end

            c = c + 1;
        end
        calculated_roots(i) = a_sol;
        fprintf('Root %5.1d: iterations = %5.1d, x_sol = %12.5f, error = %12.5f \n', i, c, a_sol, error)
    end

end

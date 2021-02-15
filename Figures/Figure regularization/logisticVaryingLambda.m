%--------------------------------------------------------------------------
% Plot the the logistic error (sigma) versus varying values of lambda.
% Lambda is the regularization parameter.
%
% INPUT FILE:
% -----------
%
%   mismatchBoolean.mat:
%   Contains a data structure array with the mismatch booleans of all
%   proteins.
%
% OUTPUT FILES:
% -------------
%
%   Plots sigma versus lambda. 
%--------------------------------------------------------------------------

clear all
close all

downloadDate = '20210120';

minLambda = -4; % -4 is a min-lambda of 10^(-4)
maxLambda = -2; % -2 is a min-lambda of 10^(-2)
N = 100;        % number of lambdas in array

proposedLambda = 1e-3;

%-------------------------------START CODE---------------------------------

lambdaArray = logspace(-4,-2,N);
numLambda = length(lambdaArray);
addpath('../../Data')
fileName = ['mismatchBooleanOverTime',downloadDate,'_t3.mat'];
load(fileName)

mismatchBoolean = mismatchBooleanOverTime(end);
proteinNames = fields(mismatchBoolean);
numProteins = length(proteinNames);

alld0 = zeros(numLambda, numProteins);
allalfa = zeros(numLambda, numProteins); 
allsigma = zeros(numLambda, numProteins);

for l = 1:numLambda
    lambda = lambdaArray(l);
    [beta, sigma] = logisticRegressionProteins(mismatchBoolean, lambda);
    alld0(l,:) = - beta(1,:) ./ beta(2,:); % d0 = -b0/b1
    allalfa(l,:) = - beta(2,:);            % alfa = -b1
    allsigma(l,:) = sigma;
end

%% Plot d0, alfa and sigma versus lambda

figure() % d0 
for p = 1:numProteins
    protein = proteinNames{p};
    subplot(5,6,p)
    h = semilogx(lambdaArray, alld0(:,p), '-or', 'MarkerSize', 1.5, 'LineWidth', 0.1);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    hold on
    xline(proposedLambda)
    hold off
    ytickformat('%.2g')
    title(protein)
end
suptitle('D0 versus Lambda')
set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 11.5 7.5])
figName = fullfile('d0_vs_lambda.pdf');
saveas(gcf, figName)

figure() % alfa
for p = 1:numProteins
    protein = proteinNames{p};
    subplot(5,6,p)
    h = semilogx(lambdaArray, allalfa(:,p), '-or', 'MarkerSize', 1.5, 'LineWidth', 0.1);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    hold on
    xline(proposedLambda)
    hold off
    ytickformat('%.2g')
    title(protein)
end
suptitle('Alfa versus Lambda')
set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 11.5 7.5])
figName = fullfile('alfa_vs_lambda.pdf');
saveas(gcf, figName)

figure() % sigma
for p = 1:numProteins
    protein = proteinNames{p};
    subplot(5,6,p)
    h = semilogx(lambdaArray, allsigma(:,p), '-or', 'MarkerSize', 1.5, 'LineWidth', 0.1);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    hold on
    xline(proposedLambda)
    hold off
    ytickformat('%.2g')
    title(protein)
end
suptitle('Sigma versus Lambda')
set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 11.5 7.5])
figName = fullfile('sigma_vs_lambda.pdf');
saveas(gcf, figName)

%% Plot some logistics at different lambda values

proteinsToShow = {'NSP9', 'Spike', 'N','NS9b'};
nP = length(proteinsToShow);

% show 5 lambda values
lambdaToShow = logspace(-4, -2, 5);
nL = length(lambdaToShow);

proteinIndecesToShow = zeros(1,nP);
for p = 1:nP
    proteinIndecesToShow(p) = find(strcmp(proteinNames,proteinsToShow{p}));
end

figure()
c = 0;
for l = 1:nL
    lambda = lambdaToShow(l);
    [beta, sigma] = logisticRegressionProteins(mismatchBoolean, lambda);
    for p = 1:nP
        c = c + 1;
        index = proteinIndecesToShow(p);
        protein = proteinNames{index};
        y = mismatchBoolean.(protein)';
        X = linspace(0,length(y)-1,length(y))';

        maxy = find(y,1,'last') + 15;
        yax = y(1:maxy);
        Xax = linspace(0,length(yax)-1,length(yax))';
        xaxis = linspace(0,length(yax)-1, 10*length(yax))';

        b = beta(:,index);
        S = predict(b, xaxis);
        S(1:2) = NaN;

        subplot(nP,nL,l + (p-1)*nL)
        plot(Xax,yax,'.r','MarkerSize',10);
        hold on
        plot(xaxis, S, 'k','LineWidth',1);
        axis([0, max(Xax),0,1.1]);
        xticks([0,max(Xax)])    
        yticks([0,1])
        
        if p == 1
            title(['Lambda = ',num2str(lambda)])
        end
        
        if l == 1
            ylabel(protein)
        end
    end
end

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 11.5 6.5])
figName = fullfile('logisticExamples.pdf');
saveas(gcf, figName)

%% Functions

function [beta, sigma] = logisticRegressionProteins(mismatchBoolean, lambda)

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

function g = sigmoid(z)
    g = 1 ./ (1 + exp(-z));
end

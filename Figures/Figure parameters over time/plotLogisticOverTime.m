%--------------------------------------------------------------------------
% Plot the parameters (d0 and alfa) of the logistic regression over time.
%
% INPUT FILE:
% -----------

%   mismatchBooleanOverTime<date>.mat:
%   Contains a data structure array with the mismatch booleans at 
%   different timepoints.

% OUTPUT FILES:
% -------------

%   Plots of d0 and alfa versus 
%--------------------------------------------------------------------------

clear all
close all

threshold = 5;
downloadDate = '20210120';
fileName = ['mismatchBooleanOverTime', downloadDate, '_t', num2str(threshold), '.mat'];
load(fileName)

pNames = fields(mismatchBooleanOverTime(1));
numProteins = length(pNames);
numTimeBins = length(timePoints) - 1;
timeCenters = mean( [timePoints(2:end); timePoints(1:end-1)] );

d0 = NaN(numTimeBins, numProteins);
alfa = NaN(numTimeBins, numProteins);
numSeqs = zeros(1, numTimeBins);

for i = 1:numTimeBins
    t2 = timePoints(i + 1);
    numSeqs(i) = sum(dates < t2);
    mismatchBoolean = mismatchBooleanOverTime(i);
    [beta, sigma] = logisticRegressionProteins(mismatchBoolean);
    d0(i,:) = - beta(1,:) ./ beta(2,:); % d0 = -b0/b1
    alfa(i,:) = - beta(2,:);            % alfa = -b1
end

%% Plot d0 and alfa over time

minD0 = min(d0(:));
maxD0 = max(d0(:));

minAlfa = min(alfa(:));
maxAlfa = max(alfa(:));

figure()
for p = 1:numProteins
    protein = pNames{p};
    subplot(5,6,p)
    h = plot(timeCenters, d0(:,p), '-or', 'MarkerSize', 1.5, 'LineWidth', 0.1);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    ylim([minD0 - 1, maxD0 + 1]);
    title(protein)
end
suptitle('d0 over time')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 11.5 7.5])
figName = fullfile('d0_over_time.pdf');
saveas(gcf, figName)

figure()
for p = 1:numProteins
    protein = pNames{p};
    subplot(5,6,p)
    h = plot(timeCenters, alfa(:,p), '-or', 'MarkerSize', 1.5, 'LineWidth', 0.1);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    ylim([minAlfa - 0.2, maxAlfa + 0.2]);
    title(protein)
end
suptitle('alfa over time')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 11.5 7.5])
figName = fullfile('alfa_over_time.pdf');
saveas(gcf, figName)

%% Plot d0 and alfa vs num sequences
figure()
for p = 1:numProteins
    protein = pNames{p};
    subplot(5,6,p)
    h = plot(numSeqs, d0(:,p), '-or', 'MarkerSize', 1.5, 'LineWidth', 0.1);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    ylim([minD0 - 1, maxD0 + 1]);
    title(protein)
end
suptitle('d0 vs nuber of sequences')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 11.5 7.5])
figName = fullfile('d0_vs_numseqs.pdf');
saveas(gcf, figName)

figure()
for p = 1:numProteins
    protein = pNames{p};
    subplot(5,6,p)
    h = plot(numSeqs, alfa(:,p), '-or', 'MarkerSize', 1.5, 'LineWidth', 0.1);
    set(h, 'MarkerFaceColor', get(h,'Color'));
    ylim([minAlfa - 0.2, maxAlfa + 0.2]);
    title(protein)
end
suptitle('alfa vs nuber of sequences')

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Units','inches','Position',[1 1 11.5 7.5])
figName = fullfile('alfa_vs_numseqs.pdf');
saveas(gcf, figName)

%%
function [beta, sigma] = logisticRegressionProteins(mismatchBoolean)

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



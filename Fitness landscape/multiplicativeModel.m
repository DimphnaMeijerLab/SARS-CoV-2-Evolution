clear all
close all

downloadDate = '20210120';
addpath('../Simulation_code')
addpath('../Data')

fitnessFunction = @truncation_model_fitness;              % d0 = -b0/b1

%% Do logistic regression

fileName = ['mismatchBooleanOverTime',downloadDate,'_t3.mat'];
load(fileName)

mismatchBoolean = mismatchBooleanOverTime(end);
proteinNames = fields(mismatchBoolean);
numProteins = length(proteinNames);

lambda = 0;
[beta, sigma] = logisticRegressionProteins(mismatchBoolean, lambda);

%% Calculate multiplicative fitness and plot

figure()
for p = 1:length(proteinNames)
    protein = proteinNames{p};

    y = mismatchBoolean.(protein)';
    N = length(y);
    X = transpose(0 : N-1);

    % Make plot
    maxy = find(y,1,'last') + 15;
    yax = y(1:maxy);
    Xax = linspace(0,length(yax)-1,length(yax))';
    xaxis = linspace(0,length(yax)-1, 10*length(yax))';
    S = fitnessFunction(beta(:,p), xaxis);
    %S(1:2) = NaN;

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

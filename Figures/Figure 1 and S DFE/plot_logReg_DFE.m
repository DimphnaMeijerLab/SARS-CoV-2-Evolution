% close all
% clear all

% Created 26 AUG 2020.
% This code plots:
%   1) Mutation number booleans and fitted logistic regressions;
%   2) Analysis of d0, alfa and sigma for each protein (how they correlate
%      and change over time);
%   3) Distribution of submissions on GISAID database;
%   4) Distribution of fitness effects.

%% Load alignment data & do logistic regression
addpath('../../Data');
addpath('../../Fitness landscape');
addpath('../../Simulation_code');

downloadDate = '20210120';
fitnessModel = 'multiplicative';
sigma = zeros(1,length(pNames)) + 0.1;

%% Define fitness function and do logistic regression

if strcmp(fitnessModel, 'sigmoid')
    fitnessFunction = @predict;
elseif strcmp(fitnessModel, 'multiplicative')
    fitnessFunction = @multiplicative_fitness;
elseif strcmp(fitnessModel, 'eigen')
    fitnessFunction = @eigen_model_fitness;
elseif strcmp(fitnessModel, 'truncation')
    fitnessFunction = @truncation_model_fitness;
else
    error('Invalid fitness model.')
end

fileName = ['mismatchBooleanOverTime',downloadDate,'_t3.mat'];
load(fileName)

mismatchBoolean = mismatchBooleanOverTime(end);
[beta, ~] = logisticRegressionProteins(mismatchBoolean, lambda);

alfa = - beta(2,:);
d0 = beta(1,:) ./ alfa;

%% Logistic regression for envelope protein

E = 16;                         % specify protein
protein = proteinNames{E};
y = mismatchBoolean.(protein)';
maxy = find(y,1,'last') + 15;   % plot up to 15 AAs after the last observed hamming distance
y = y(1:maxy);
X = linspace(0,length(y)-1,length(y))';
b = beta(:,E);
xaxis = linspace(0,length(y)-1, 10*length(y))';

close all
figure()
plot(X,y,'.r','MarkerSize',10);
hold on
plot(xaxis, fitnessFunction(b, xaxis), 'k','LineWidth',1);

xlabel('Hamming distance, d')
ylabel('Expected fitness \newline contribution, <w_{Spike}>')

yticks([0,1])
ylim([-0.05, 1.05]);

plot([d0(E) d0(E)],[-10 0.5],'k:','LineWidth',2)
hold on
plot([0 d0(E)],[0.5 0.5],'k:','LineWidth',2)
yticks(0:0.5:1)
yticklabels([0 0.5 1])

leg = legend('Spike','Logistic fit','Location','east');
leg.ItemTokenSize = [5, 1];

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 2.165 1.74])
saveas(gcf,'Figures/Spike.pdf')

%% Fitness contribution of one protein

figure()
r_expected = 0.3;
x = 0:0.01:100;
pd = makedist('Normal', r_expected, sigma(E));
t = truncate(pd,0,Inf);
rdisttrib = pdf(t,x)/sum(pdf(t,x));
plot(x,rdisttrib,'k','LineWidth',1)
hold on
area(x,rdisttrib,'FaceAlpha',0.2,'FaceColor',[0 0 0], 'EdgeColor',[0 0 0])
% hold on
% plot([r_expected r_expected],[0 max(rdisttrib)],'k','LineWidth',1)

xlim([0 1.5])
yticklabels([])
xlabel('Fitness contribution, w_{Spike}')
ylabel('Probability')

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[4 1 2.165 1.74])
saveas(gcf,'Figures/Distribution_Spike.pdf')

%% Logistic regression for NSP9
E = 9;
protein = proteinNames{E};
y = mismatchBoolean.(protein)';
maxy = find(y,1,'last') + 15;
y = y(1:maxy);
X = linspace(0,length(y)-1,length(y))';
b = beta(:,E);
xaxis = linspace(0,length(y)-1, 10*length(y))';

figure()
plot(X,y,'.r','MarkerSize',10);
hold on
plot(xaxis, fitnessFunction(b, xaxis), 'k','LineWidth',1);
xlabel('Hamming distance, d')
ylabel('Expected fitness \newline contribution, <w_{NSP9}>')
yticks([0,1])
ylim([-0.05, 1.05]);
yticks(0:0.5:1)
yticklabels([0 0.5 1])
leg = legend('NSP9','Logistic fit','Location','east');
leg.ItemTokenSize = [5, 1];

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[5 1 2.165 1.74])
saveas(gcf,'Figures/NSP9.pdf')

%% Plot all logistic regressions in one figure
figure()
proteinNamesplot = proteinNames;
proteinNamesplot{16} = 'S';

for E = 1:length(proteinNames)
    protein = proteinNames{E};
    y = mismatchBoolean.(protein)';
    X = linspace(0,length(y)-1,length(y))';
    %a = glmfit(X,y,'binomial','link','logit');
    
    maxy = find(y,1,'last') + 15;
    yax = y(1:maxy);
    Xax = linspace(0,length(yax)-1,length(yax))';
    xaxis = linspace(0,length(yax)-1, 10*length(yax))';
    b = beta(:,E);
    S = fitnessFunction(b, xaxis);
    S(1:2)=NaN;
 
    subplot(5,6,E)
    plot(Xax,yax,'.r','MarkerSize',10);
    hold on
    plot(xaxis, S, 'k','LineWidth',1);
    axis([0, max(Xax),0,1.1]);
    xticks([0,max(Xax)])    

    yticks([0,1])
    title(proteinNamesplot{E})
end

xlabel('Hamming Distance, d')
%ylabel('1(d\in GISDAID)','Color','r')
set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[5 3 3*2.165 2*1.74])
saveas(gcf,'Figures/Logistics.pdf')

%% Plot fitness as heatmap

Nmax = 40; % maximum hamming distance to plot
x = 0:Nmax;
rAxis = linspace(0, 2.5, 100);
% figure()
% for i = 1:length(pNames)
%     
%     beta = [beta0(i), beta1(i)]';
%     mu_array = predict(beta, x');
% 
%     map = zeros(length(rAxis), Nmax);
%     for j = 1:length(mu_array)
%         pd = makedist('Normal', mu_array(j), sigma(i));
%         t = truncate(pd, 0, inf);
%         map(:, j) = pdf(t, rAxis);
%     end
%     
%     subplot(6,4,i)
%     imagesc('XData', x, 'YData', rAxis, 'CData', map / max(map(:)))
%     ylim([0,1.5])
%     yticks(0:.5:1.5)
%     xlim([0, Nmax])
%     colormap summer
%     hold on
%     xaxis = linspace(0, Nmax, 10*Nmax)';
%     plot(xaxis, predict(beta, xaxis), '-r', 'LineWidth', 2)
%     hold off
%     title(pNames{i})
% end
% 
% set(gcf,'Color','w','Position',[450 500 400 500],'Units','inches')
% %saveas(gcf,'Figures/replRate.pdf')


%% Histogram of points per time
%{
figure()
histogram(dates, time,'FaceColor','k','FaceAlpha',1)
ylabel('Frequency')

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[8 1 2*2.165 1.74])
saveas(gcf,'Figures/timePoints.pdf')

%% Store all alfa's and d0's in an array

allAlfa = zeros(length(proteinNames), length(time)-1);
allD0 = zeros(length(proteinNames), length(time)-1);
for t = 1:length(time)-1
    allAlfa(:, t) = - mismatchesPerTime(t).beta1';
    allD0(:, t) = mismatchesPerTime(t).beta0' ./ allAlfa(:, t);
end

leg = cell(1, length(proteinNames));
for p = 1:length(proteinNames)
    leg{p} = proteinNames{p};
end

%% Normalize max = 1
for p = 1:length(proteinNames)
    allAlfa_max1(p, :) = allAlfa(p, :) / max( allAlfa(p, :) );
    allD0_max1(p, :) = allD0(p, :) / max( allD0(p, :) );
    leg{p} = proteinNames{p};
end
% figure()
% subplot(2,1,1)
% plot(time(2:end), allD0_max1, '-o')
% ylabel('d_0')
% legend(leg)
% 
% subplot(2,1,2)
% plot(time(2:end), allAlfa_max1, '-o')
% ylabel('steppness \alfa','Interpreter', 'Latex')
% legend(leg)
% 
% set(gca, 'FontSize',9, 'LineWidth',1)
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'Color','w','Position',[250 250 800 500],'Units','inches')
% saveas(gcf,'Figures/alfa_d0_normalized_max_1.pdf')

% Mean d0
meanD0 = mean(allD0_max1);
stdD0 = 1.96 * std(allD0_max1) / sqrt(length(proteinNames));
% figure()
% errorbar(datenum(time(2:end)), meanD0, stdD0)
% ylabel('mean d0')
% xlabel('time')

%% Plot beta1 and beta2 versus protein length

colors = {'r','g','k','c','m','k','w'};
markers = {'o','+','*','x','s','d'};

c=0;
lengthsArray = [];
beta1Array = [];
maxL = 0;
leg = {};
% figure()
% for p=1:length(proteinNames)
%     proteinName = proteinNames{p};
%     L = pInfo.(proteinName).La; 
%     if L > maxL
%         maxL = L;
%     end
%     
%     Ncolor = max(1,mod(p, length(colors)));
%     Nmarker = ceil(p / length(colors));
%     colorMarker = [colors{Ncolor}, markers{Nmarker}];
%     subplot(2,1,1)
%     plot(L, alfa(p), colorMarker)
%     hold on
%     
%     c = c + 1;
%     lengthsArray(c) = L;
%     colorMarkers{c} = [colors{Ncolor}, markers{Nmarker}];
% 
%     subplot(2,1,2)
%     plot(L, d0(p), colorMarkers{c})
%     leg{c} = proteinName;
%     hold on
% end
% subplot(2,1,1)
% legend(proteinNames)
% title('alfa, steepness, versus protein length')
% xlabel('Protein length (number of amino acids)')
% ylabel('\alfa')
% 
% subplot(2,1,2)
% legend(leg)
% title('\d0, inflection point, versus protein length')
% xlabel('Protein length (number of amino acids)')
% ylabel('d0')
% 
% suptitle('alfa and d0 versus')
% 
% saveas(gcf,'Figures/BetaAnalysis.fig')

%% d0, alfa and sigma per protein

[d0_sorted, ind] = sort(d0);
alfa_sorted = alfa(ind);
sigma_sorted = sigma(ind);

pNames_sorted = pNames(ind)' ;

AllParameters_sorted = table(pNames_sorted , d0_sorted', alfa_sorted',sigma_sorted')

figure()
subplot(3,1,1)
bar(d0_sorted,'EdgeColor',[0,0,0],'FaceAlpha',0.5,'FaceColor',[0 0 0])
ylim([-5 max(d0_sorted)*1.1])
ylabel('$d_0$','Interpreter', 'latex')
xticks(1:24)
xticklabels(proteinNames(ind))
xtickangle(45)

subplot(3,1,2)
bar(alfa_sorted,'EdgeColor',[0,0,0],'FaceAlpha',0.5,'FaceColor',[0 0 0])
ylim([0 max(alfa_sorted)*1.1])
ylabel('\alpha')
xticks(1:24)
xticklabels(proteinNames(ind))
xtickangle(45)

subplot(3,1,3)
bar(sigma_sorted,'EdgeColor',[0,0,0],'FaceAlpha',0.5,'FaceColor',[0 0 0])
ylim([0 max(sigma_sorted)*1.1])
ylabel('\sigma')
xticks(1:24)
xticklabels(proteinNames(ind))
xtickangle(45)
saveas(gcf,'Figures/LogisticsData.pdf')

%% Time course of GISAID database
for i=1:length(time)
   NumSeqs(i) = length(find(datenum(dates)<datenum(time(i)))); 
end

figure()
plot(time,NumSeqs,'k','LineWidth',1)
xlabel('Date')
ylabel('Sequence count')
yticks((0:2:10)*10^4)

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[8 1 3*2.165 1.74])
saveas(gcf,'Figures/NumSeqsvsTime.pdf')
%}

%% Distribution of fitness effects

% Load simulated data
load(['Data/simulatedDFE_',fitnessModel,'.mat'])
r0 = 1.5;

colorbar = jet;

figure()
h = histogram(r(1,:)/r0,'FaceColor',colorbar(1,:),'Normalization','probability');
%hold on
X = h.BinEdges((1:h.NumBins)+1)-h.BinWidth/2;
Y = h.Values/sum(h.Values);
plot(X,Y,'Color',colorbar(1,:))
hold on
area(X,Y,'FaceAlpha',0.2,'FaceColor',colorbar(1,:))

xlabel('Fitness, W')
ylabel('Probability')
xlim()
xticks(0:0.25:1.25)
ylim([0 0.2])
yticks(0:0.05:0.2)
set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 2.165 1.74])
saveas(gcf,'Figures/DFE1.pdf')

fig = figure();
%subplot(2,1,2)
v = violinplot(r' / r0);
% Specify color of violins

for d = 1:maxD
    v(d).ShowData = false;
    v(d).ViolinColor = colorbar(5*d,:); %[0.5 0.5 0.5];
    v(d).BoxColor = [0 0 0];
    v(d).ShowMean = false;
end

yline(1);
%title('Normal distribution')
xlabel('Hamming distance, d')
ylabel('Fitness, W')

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 1.8*2.165 1.74])
fig.Renderer = 'Painters';
saveas(gcf,'Figures/ViolinplotFitnessLandscape.pdf')

%% Table
alfa = - beta(2,:);
d0 = - beta(1,:) ./ beta(2,:);

C = cell(length(pNames), 4);
for p = 1:length(pNames)
    protein = pNames{p};
    C{p,1} = protein;
    C{p,2} = length(mismatchBoolean.(protein));
    C{p,3} = d0(p);
    C{p,4} = alfa(p);
    C{p,5} = sigma(p);
    C{p,6} = nSeqPerProtein.(protein);
end
T = cell2table(C,'VariableNames',{'name','length','d0','alfa','sigma','nSeqs'});
writetable(T,'table1.csv')


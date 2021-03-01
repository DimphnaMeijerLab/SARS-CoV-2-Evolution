clear all
close all

%% Collect all data
folder = 'Output';
files = dir([folder, '/Gillespie*']);
Nsim = length(files);
alldata = struct();

for vb = 1:Nsim
    fname = [folder,'/', files(vb).name];
    load(fname);
    
    alldata(vb).('time') = data.t;
    alldata(vb).('ntot') = data.ntot;
    alldata(vb).('nAA') = data.nAA;
    alldata(vb).('U') = data.U;
    alldata(vb).('I') = data.I_sum;
    alldata(vb).('V') = data.V_sum;
    alldata(vb).('maxD') = data.maxD;
    alldata(vb).('maxR') = data.maxR;
    alldata(vb).('statY') = data.statY;
    alldata(vb).('relativeY') = data.relativeY;
    alldata(vb).('statR') = data.statR;
    alldata(vb).('statD') = data.statD;
    alldata(vb).('diversity') = data.diversity;
end

% Collect all time points in one array:
alltime = vertcat(alldata(:).time);
[sortedTime, ind] = sort(alltime);
Npoints = Nsim;

% Define colors:
shadyred = [250,189,189];
shadyred = shadyred / max(shadyred);
red = {shadyred, 'r'};

shadyblue = [196, 201, 239];
shadyblue = shadyblue / max(shadyblue);
blue = {shadyblue, 'b'};

lightgrey = [224, 224, 224];
lightgrey = lightgrey / 255;
grey = {lightgrey, 'k'};

orange = [204 102 0];
orange = orange / max(orange);

%% figure 1: time course of infection (mean and SD)
close all
% Calculate mean and error of V, U and I:
allV = vertcat(alldata(:).V);
sortedV = allV(ind);
[timeAxis_V, V_mean, V_error] = timeAverage(sortedTime, sortedV, Npoints);

allI = vertcat(alldata(:).I);
sortedI = allI(ind);
[timeAxis_I, I_mean, I_error] = timeAverage(sortedTime, sortedI, Npoints);

allU = vertcat(alldata(:).U);
sortedU = allU(ind);
[timeAxis_U, U_mean, U_error] = timeAverage(sortedTime, sortedU, Npoints);

% Plot mean with shady error
figure(1)
h = gobjects(3,1);
h(1) = plotShadyError(timeAxis_U, U_mean, U_error, blue,'log');
h(2) = plotShadyError(timeAxis_I, I_mean, I_error, red,'log');
h(3) = plotShadyError(timeAxis_V, V_mean, V_error, grey,'log');
set(gca, 'YScale', 'log' )
leg = legend(h, 'Non-infected cells', 'Infected cells', 'Free viral particles');
leg.ItemTokenSize = [5, 1];

xlabel('Time (days p.i.)')
ylabel('Number')
yticks(10.^[0:6])
set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 2*2.165 1.74])
saveas(gcf,'Figures/TimeCourseGraph.pdf')

%% figure 2: time course of infection (all simulations)

% figure(2)
% 
% for vb = 1:Nsim
%     t = alldata(vb).('time');
%     U = alldata(vb).('U');
%     V = alldata(vb).('V');
%     I = alldata(vb).('I');
%     
%     semilogy(t, U,'b', 'LineWidth',2)
%     hold on
%     semilogy(t, V,'r', 'LineWidth',2)
%     hold on
%     semilogy(t, I,'Color',orange,'LineWidth',2)
%     
%     hold on
% end
% legend({'Non-infected cells','Free viral particles','Infected cells'},'FontSize',9)
% leg.ItemTokenSize = [15, 1];
% ylabel('Number')
% xlabel('Time (days p.i.)')
% set(gca,'LineWidth',1)
% 
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'Color','w','Position',[200 200 900 550],'Units','inches')

%% Figure 3: Number of strains

figure(3)
allntot = vertcat(alldata(:).ntot);
allnAA = vertcat(alldata(:).nAA);

sortedntot = allntot(ind);
sortednAA = allnAA(ind);
[timeAxis_ntot, ntot_mean, ntot_error] = timeAverage(sortedTime, sortedntot, Npoints);
[timeAxis_nAA, nAA_mean, nAA_error] = timeAverage(sortedTime, sortednAA, Npoints);

h = gobjects(2,1);
h(1) = plotShadyError(timeAxis_ntot, ntot_mean, ntot_error, grey,'log');
h(2) = plotShadyError(timeAxis_nAA, nAA_mean, nAA_error, blue,'log');

leg = legend(h, 'Distinct nucleotide sequences', 'Distinct amino acid sequences')
leg.ItemTokenSize = [5, 1];

ylabel('Number')
xlabel('Time (days p.i.)')
yticks(10.^[0:5])

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 2*2.165 1.74])
saveas(gcf,'Figures/TimeCourseNumberofSequences.pdf')

%% Figure 4: maximal reproduction rate over time 

allmaxR = horzcat(alldata(:).maxR)';
sortedMaxR = allmaxR(ind);

Npoints = Nsim;
[timeAxis, maxR_mean, maxR_error] = timeAverage(sortedTime, sortedMaxR, Npoints);

figure(4)

p = plotShadyError(timeAxis, maxR_mean / params.r0, maxR_error, red,'normal');
hold on
yline(mean(horzcat(alldata(:).statR))/params.r0,'k--','Alpha',1,'LineWidth',2)
ylabel('Maximum fitness')
xlabel('Time (days p.i.)')
ylim([0.95 1.5])
yticks(1:0.1:1.5)

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 2.165 1.74])
saveas(gcf,'Figures/MaxRepRate.pdf')


%% Figure 5: abundance of hamming distances

% Find the maximal hamming distance observed in all simulations (maxY):
maxY = 0;
for vb = 1:Nsim
    lengthY = length(alldata(vb).statY);
    if maxY < lengthY
        maxY = lengthY;
    end
end

% Collect all stationary Y values
statY = zeros(Nsim, maxY);
for vb = 1:Nsim
    if sum(isnan(alldata(vb).statY)) == 0
        statY_vb = alldata(vb).statY;
        statY_vb = [statY_vb, zeros(1, maxY - length(statY_vb))]; % add zeros to make the right length;
        statY(vb, :) = statY_vb;
    end
end

% Plot histogram
figure(5)
barcolor = [.2 .6 .5];
maxY = find(sum(statY)>0, 1, 'last');
statY = statY(:, 1:maxY);
x = linspace(0, maxY-1, maxY);

meanY = mean(statY);
errorY = std(statY);
er = errorbar(x,meanY, errorY);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on
bar(x, meanY,'FaceColor', barcolor)
set(gca,'yscale','log')
xticks([0:(maxY-1)])

% for i=1:length(x)
%     text(x(i),meanY(i)*2,strcat(num2str(meanY(i),2),'%'),'HorizontalAlignment','center')
% end


%title('Relative abundance of hamming distances')
ylim([10^(-12) 10])
yticks(10.^(-12:2:1))

xlabel('Hamming distance')
ylabel('Stationary fraction')

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 2.165 1.74])
saveas(gcf,'Figures/Ndistribution.pdf')


%% Figure 6: distribution of max values
figure(6)
maxR_all = zeros(Nsim,1); 
maxD_all = zeros(Nsim,1);
statR_all = zeros(Nsim,1);
statD_all = zeros(Nsim,1);

for i = 1:Nsim
    maxR_all(i) = max(alldata(i).maxR);
    maxD_all(i) = max(alldata(i).maxD);
    statR_all(i) = alldata(i).statR;
    statD_all(i) = alldata(i).statD;
end
subplot(2,1,1)
histogram(maxR_all)
xlabel('Max replication rate during the infection')
ylabel(['Frequency (out of ', num2str(Nsim),' simulations)'])

subplot(2,1,2)
histogram(maxD_all)
xlabel('Max hamming distance during the infection')
ylabel(['Frequency (out of ', num2str(Nsim),' simulations)'])

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Position',[200 200 900 550],'Units','inches')
%saveas(gcf,'Figures/maxRD_dist.pdf')

%% Figure 7: distribution of stationary values
figure(7)
subplot(2,1,1)
histogram(statR_all)
xlabel('Stationary replication rate during the infection')
ylabel(['Frequency (out of ', num2str(Nsim),' simulations)'])

subplot(2,1,2)
histogram(statD_all)
xlabel('Stationary hamming distance during the infection')
ylabel(['Frequency (out of ', num2str(Nsim),' simulations)'])

set(gcf,'PaperOrientation','landscape');
set(gcf,'Color','w','Position',[200 200 900 550],'Units','inches')
%saveas(gcf,'Figures/statRD_dist.pdf')

%% Figure 8: heatmap

allRelY = zeros(length(alltime), maxY);
c = 1;
for i = 1:Nsim
    relativeY = alldata(i).relativeY;
    [n, m] = size(relativeY);
    relativeY = [relativeY, zeros(n, maxY - m)]; % add zeros to make the right length;]
    allRelY(c:c+n-1, :) = relativeY;
    c = c + n;
end

sortedAllRelY = allRelY(ind, :);
[timeAxis_relY, relY_means, relY_error] = timeAverage(sortedTime, sortedAllRelY, Npoints);

figure(8)

%subplot(2,1,1)
HeatMap = imagesc('XData', timeAxis_relY, 'CData', relY_means');
set(gca, 'ColorScale', 'log')
colormap(jet)
z = colorbar;

%alpha(HeatMap,0.5);

xlim([timeAxis_relY(1), timeAxis_relY(end)])
maxY = size(relY_means,2);
yticks(1:maxY)
ylabels = cell(1, maxY);
for i = 1:maxY
    ylabels{i} = num2str(i - 1);
end
yticklabels(ylabels)
ylabel('Hamming distance')
xlabel('Time (days p.i.)')
z.Label.String = 'Relative abundance';
z.Ticks = 10.^(-5:1);
z.Direction = 'reverse'; 

% % Plot error
% subplot(2,1,2)
% imagesc('XData', timeAxis_relY, 'CData', relY_error');
% set(gca, 'ColorScale', 'log')
% colormap(summer)
% z = colorbar;
% 
% xlim([timeAxis_relY(1), timeAxis_relY(end)])
% maxY = size(relY_means,2);
% yticks(1:maxY)
% yticklabels(ylabels)
% ylabel('Hamming distance')
% xlabel('time (days p.i.)')
% z.Label.String = 'Error';

set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 2*2.165 1.74])
saveas(gcf,'Figures/HeatMap.pdf')

%% Figure 9 + 10: Relative abundance of Hamming distances over time

[T, maxD] = size(relY_means);
figure(9)
b = bar(timeAxis_relY, fliplr(relY_means), 1, 'stacked');
set(gca, 'YScale', 'log')
barPositions = zeros(length(b)+1, length(timeAxis_relY));
timeAxis = [timeAxis_relY', fliplr(timeAxis_relY')];

colorbar = jet;

for d = 1:length(b)
    barPositions(d+1, :) = b(d).YData;
    lgd{d} = ['d = ', num2str(d - 1)];
end

barPositions(barPositions < 1e-9) = 1e-9;
barPositions = flipud(barPositions);

figure(10)
for d = 1:length(b)
    curve = [barPositions(d+1, :), fliplr(barPositions(d, :))];
    fill(timeAxis, curve, colorbar(7*d,:), 'LineStyle', 'none')
    hold on
end
set(gca, 'YScale', 'log')
legend(lgd,'Location','eastoutside')
xlabel('Time (days p.i.)')
ylabel('Relative abundance')
yticks(10.^[-10:2:0])
xlim([0 20])
set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 2*2.165 1.74])
saveas(gcf,'Figures/RelAbTime.pdf')

%% Figure 11: diversity index

allDiversity = horzcat(alldata(:).diversity)';
sortedDiversity = allDiversity(ind);
nans = isnan(sortedDiversity);

sortedDiversity = sortedDiversity(~nans);
sortedTime_div = sortedTime(~nans);

[timeAxis_div, div_mean, div_error] = timeAverage(sortedTime_div, sortedDiversity, Npoints);

figure(11)
p = plotShadyError(timeAxis_div, div_mean, div_error, grey, 'normal');
xlabel('Time (days p.i.)')
ylabel('Diversity') %Gini-Simpson diversity index
ylim([-0 0.65])
yticks(0:0.1:0.6)
xlim([0 20])
set(gca, 'FontSize',9, 'LineWidth',1)
set(gcf,'Color','w','Units','inches','Position',[1 1 2.165 1.74])
saveas(gcf,'Figures/Diversity.pdf')


%% Functions

function [timeAxis, means, error] = timeAverage(t, y, Npoints)
% This function takes a time average of y.
% y is an nxm array, with n the number of time points.
% time is a column vector with length n.
% Npoints is the number of points used to take an average. It can be set to
% the number of simulations.

    [n,m] = size(y);
  
    time = [t(1:Npoints:end); max(t)];
    timeAxis = (time(2:end) + time(1:end-1))/2;
    
    means = zeros(length(timeAxis), m);
    error = zeros(length(timeAxis), m);

    for i = 1:length(timeAxis)
        inBin = (t >= time(i) & t < time(i+1));
        means(i, :) = nanmean(y(inBin, :));
        error(i, :) = nanstd(y(inBin, :));
    end
end

function p = plotShadyError(x, means, error, color, yscale)

    if strcmp(yscale, 'log')
        curve1 = max(means + error, 1);
        curve2 = max(means - error, 1);
        xaxis = [x; flipud(x)];
        inBetween = [curve1; flipud(curve2)];

        fill(xaxis, inBetween, color{1});
        hold on;

        means(means < 1) = 1;
        p = plot(x, means, color{2}, 'LineWidth', 1);
        set(gca, 'YScale', 'log' )

    else
        curve1 = means + error;
        curve2 = means - error;
        xaxis = [x; flipud(x)];
        inBetween = [curve1; flipud(curve2)];

        fill(xaxis, inBetween, color{1});
        hold on;

        p = plot(x, means, color{2}, 'LineWidth', 1);
        set(gca, 'YScale', 'default')
    end
end

function [t, x] = solveDE(tRange, initialCond, params)
    
    a = params(1);
    b = params(2);
    c = params(3);
    r0 = params(4);
    
    dfdt = @(t,x) [r0 * x(2) - b * x(1) - a * x(3) * x(1); ... 
                   a * x(3) * x(1) - c * x(2); ... 
                   - a * x(3) * x(1)];
                       
    [t,x] = ode45(dfdt, tRange, initialCond);
end
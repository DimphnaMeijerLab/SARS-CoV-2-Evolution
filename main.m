close all
addpath('./Simulation_code');
addpath('./Data');

% Parameter initialisation
nr = 1;             % simulation number
U0 = 1e4;           % initial number of uninfected cells
mu_array = 1e-6;    % mutation rate(s) to test

% Run simulation
[data, params] = gillespie(nr, U0, mu_array);

% Plot outcome
t = data.t;
U = data.U;
I = data.I_sum;
V = data.V_sum;

semilogy(t, U, 'k', t, I, 'b', t, V, 'r')
legend('Uninfected cells', 'Infected cells', 'Free viral particles')
xlabel('Time (days p.i.)')
ylabel('Count')

%% SA
PRCC_mat, pval_mat = sensitivity_analysis(10, 2);

% plot
data_fields_plot = {'Peak time', 'Peak load', 'End time', 'stat(r)', 'max(r)', 'stat(diversity)', 'stat(d)'};
fig1 = figure();
movegui(fig1,'west');

PRCC_mat_plot = PRCC_mat;
PRCC_mat_plot(~(pval_mat<0.05)) = nan;
hm = heatmap(round(PRCC_mat_plot,2),'ColorMap', redblue(),'MissingDataLabel','n.s.','MissingDataColor',[0.4 0.4 0.4]);
hm.XDisplayLabels = params_fields;
%hm.XLabel = 'Input parameters';
hm.YDisplayLabels = data_fields_plot;
%hm.YLabel = 'Output responses';
%hm.Title = 'PRCC matrix';
caxis(hm,[-1 1])


set(gca, 'FontSize',9)
set(gcf,'Color','w','Units','inches','Position',[1 1 2*2.165 1*1.74]) 
% saveas(fig1, ['Figures', filesep, 'PRCC_heatmap.pdf'])

%% 
function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b];
end
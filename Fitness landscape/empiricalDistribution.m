clear all
close all

% Probabilities of deleterious versus beneficial mutations
p_deleterious = 0.958;
p_beneficial = 0.042;

% Number of simulations
N = 1e6;
N_deleterious = round(p_deleterious * N);
N_beneficial = round(p_beneficial * N);

% Make probability distributions and truncate them
pd_deleterious = makedist('Lognormal','mu',log(0.092),'sigma',1.206);
pd_deleterious = truncate(pd_deleterious,0,1);
pd_beneficial = makedist('Exponential','mu',1);
pd_beneficial = truncate(pd_beneficial,1,inf);

% Draw fitnesses from these distributions
f_deleterious = random(pd_deleterious, 1, N_deleterious);
f_beneficial = random(pd_beneficial, 1, N_beneficial);
f = [f_deleterious,f_beneficial];

% Make histogram
figure()
subplot(2,2,1)
h = histogram(f_deleterious, 'normalization', 'probability', 'FaceColor', 'k');
title('Deleterious mutations')
subplot(2,2,2)
histogram(f_beneficial, 'normalization', 'probability')
title('Beneficial mutations')
subplot(2,2,3:4)
histogram(f,'normalization', 'probability')
title('All mutations')
%set(gca, 'YScale', 'log')

figure()
binCenters = (h.BinEdges(1:end-1) + h.BinEdges(2:end)) / 2;
bar(binCenters,cumsum(h.Values))
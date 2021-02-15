% This code simulates the fitness landscape.
clear all
close all
clc

%% Specification of simulation
nSim = 1e4;         % number of simulations
r0 = 1.5;           % replication rate of WT sequence
maxD = 12;          % maximal hamming distance
downloadDate = '20210120';
lambda = 1e-3;

%% Get reference sequences, beta, sigma
addpath('../../Data')
addpath('../../Simulation_code')

fileName = ['mismatchBooleanOverTime',downloadDate,'_t3.mat'];
load(fileName)
mismatchBoolean = mismatchBooleanOverTime(end);
refseq = load('refSeq.mat');
pRefSeq = refseq.pRefSeq;
pNames = refseq.pNames;
proteinLocation = refseq.proteinLocation;

[beta, sigma] = logisticRegressionProteins(mismatchBoolean, lambda);

La = length(pRefSeq);
nP = length(pNames);

%% Simulate
r = zeros(maxD, nSim);

for d = 1:maxD
    % Draw random mutation locations (uniformly distributed)
    mutLocs = randi(La, d, nSim);
    % Find out which proteins are mutated 
    mutations = zeros(d, nSim);
    for p = 1:nP
        startP = proteinLocation(p, 1);
        endP = proteinLocation(p, 2);
        mutations(mutLocs >= startP & mutLocs <= endP) = p;
    end
    % Draw replication rates 
    for i = 1:nSim
        di = zeros(nP, 1);

        for j = 1:d
            iMut = mutations(j, i);
            di(iMut) = di(iMut) + 1;
        end

        r(d, i) = replicationRate(di, r0, sigma, beta, 'normal');
    end
end

save('Data/simulatedDFE.mat', 'r', 'nSim', 'maxD', 'r0')


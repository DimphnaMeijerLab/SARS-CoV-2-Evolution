function [PRCC_mat, pval_mat] = sensitivity_analysis(nPoints, nIter)
% Sensitivity analysis function to run parameter sets x times and average
% results and PRCC

% Latin Hypercube Sampling (LHS) intitialization
U0 = 1e4;
params_lhs = lhs_initialization(nPoints, 0.2, [4.5e-3, 0.9, 1.5, 1e-6]); % [a b r0 mu]

% Define input parameter and output response fields
params_fields = {'a', 'b', 'r0', 'mu'};
data_fields = {'V_peakTime', 'V_peak', 't_end', 'statR', 'maxmaxR', 'statDiv', 'statD'};
nField = length(data_fields);

data_cell = cell(nPoints, nIter);
params_cell = cell(nPoints, nIter);
nr = 0;
for iPoints = 1:nPoints
    for iIter = 1:nIter
        nr = nr + 1;
        [data, params] = gillespie(nr, U0,...
                                   params_lhs(iPoints, 4), ...
                                   'a', params_lhs(iPoints, 1),...
                                   'b', params_lhs(iPoints, 2),...
                                   'c', params_lhs(iPoints, 2),...
                                   'r0', params_lhs(iPoints, 3));
        data.maxmaxR = max(data.maxR);
        for iField = 1:nField
            data_filter.(data_fields{iField}) = data.(data_fields{iField}); 
        end
        data_cell{iPoints, iIter} = data_filter;
        params_cell{iPoints, iIter} = params;
    end
end
data = [data_cell{:}];
params = [params_cell{:}];
nParams = length(params);

%% Process data
% Get unique input parameters
params_mat = [[params.a]', [params.b]', [params.r0]', [params.mu]'];
[params_unique, ~, params_index] = unique(params_mat, 'rows');
nPoint = size(params_unique, 1);

% For data fields of interest mean/std/sem output responses for each unique
% parameter point
data_mat = zeros(nParams, nField);
data_mean = zeros(nPoint, nField);
data_std = zeros(nPoint, nField);
data_sem = zeros(nPoint, nField);
for iField=1:length(data_fields)
    data_mat(:,iField) = [data.(data_fields{iField})];
    data_mean(:,iField) = accumarray(params_index, data_mat(:,iField), [], @mean);
    data_std(:,iField) = accumarray(params_index, data_mat(:,iField), [], @std);
    data_sem(:,iField) = accumarray(params_index, data_mat(:,iField), [], @(x)std(x)/sqrt(length(x)));
end
          
% PRCC
[PRCC_mat, pval_mat] = partialcorri(data_mean, params_unique, 'type', 'Spearman');

%% Save
save('SA_lhs_nK_50.mat', ...
     'params_mat', 'params_index', 'params_fields', 'data_mat', 'data_fields', ... % data per simulation
     'params_unique', 'data_mean', 'data_std', 'data_sem', ... % data per point in hypercube
     'PRCC_mat', 'pval_mat') % global sensitivity analysis results
end
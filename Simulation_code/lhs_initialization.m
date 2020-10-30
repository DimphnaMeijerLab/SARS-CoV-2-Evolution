function params_lhs = lhs_initialization(nPoints, varRatio, defVal)
% Latin Hypercube Sampling (LHS) intitialization

% LHS with 1000 points for 4 parameters
nParameters = 4;
params_lhs = lhsdesign(nPoints, nParameters); % [a, b, r0, mu]

% set standard parameter values [4.5e-3, 0.9, 1.5, 1e-6]
a = defVal(1);
b = defVal(2);
r0 = defVal(3);
mu = defVal(4);

% Convert points in hypercube to input parameter sets
params_lhs = params_lhs * 2 * varRatio + (1 - varRatio); % variation around standard parameter value
params_lhs(:,1) = params_lhs(:,1) * a;
params_lhs(:,2) = params_lhs(:,2) * b;
params_lhs(:,3) = params_lhs(:,3) * r0;
params_lhs(:,4) = params_lhs(:,4) * mu;

% Save input parameter sets
fileString = sprintf('Data/params_lhs_%g_%d', nPoints, nParameters);
disp(['Saving ', fileString])
save(fileString, 'params_lhs')

end
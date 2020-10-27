% Latin Hypercube Sampling (LHS) intitialization

% LHS with 1000 points for 4 parameters
params_lhs = lhsdesign(1e3, 4); % [a b r0 mu]

% set standard parameter values
a = 4.5e-3;
b = 0.9;
r0 = 1.5;
mu = 1e-6;

% Convert points in hypercube to input parameter sets
varRatio = 0.2; % 20% variation around standard parameter value
params_lhs = params_lhs * 2 * varRatio + (1 - varRatio);
params_lhs(:,1) = params_lhs(:,1) * a;
params_lhs(:,2) = params_lhs(:,2) * b;
params_lhs(:,3) = params_lhs(:,3) * r0;
params_lhs(:,4) = params_lhs(:,4) * mu;

% Save input parameter sets
save('Data/params_lhs_1e3_4.mat', 'params_lhs')

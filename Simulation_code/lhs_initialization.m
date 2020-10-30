function params_lhs = lhs_initialization(nPoints, varRatio, defVal)
        %------------------------------------------------------------------
        % Latin Hypercube Sampling (LHS) intitialization.
        
        % INPUT PARAMETERS
        %-----------------
        
        % nPoints: int
        % number of points in the 4D hypercube.
        
        % varRatio: double
        % Specifies the variation of the parameter values around the 
        % default values as a ratio of the default values. 
        
        % defVal: array of 4 doubles
        % Default values around which values are varied in the hypercube.
        
        % OUTPUT PARAMETERS
        %-----------------
        
        % params_lhs: nPoints x 4 double array containing coordinates of
        % points in the 4D hypercube, which correspond to parameter values
        % of the infection rate (a), clearance/death rate (b/c), the
        % reference replication rate (r0) and the mutation rate (mu).
        
        %------------------------------------------------------------------

% LHS with nPoints parameter sets for 4 parameters
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
addpath('./Simulation_code');
addpath('./Data');

% Parameter initialisation
nr = '1';           % simulation number
U0 = 1e4;           % initial number of uninfected cells
mu_array = 1e-6;    % mutation rate(s) to test

% Run simulation
[data, params] = gillespie(nr, U0, mu_array);

% Plot outcome
data_collect = data.data_collect;
t = data_collect(:,1);
U = data_collect(:,4);
I = data_collect(:,5);
V = data_collect(:,6);

semilogy(t, U, 'k', t, I, 'b', t, V, 'r')
legend('Uninfected cells', 'Infected cells', 'Free viral particles')
xlabel('Time (days p.i.)')
ylabel('Count')
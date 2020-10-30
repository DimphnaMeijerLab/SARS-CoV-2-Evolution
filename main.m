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

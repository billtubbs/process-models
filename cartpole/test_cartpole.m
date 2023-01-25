% test the following files:
%
%  - cartpole_xkp1.m
%  - cartpole_yk.m
%
% TODO: Add tests for these files:
%  - pend_StateFcn
%  - pend_MeasurementFcn.m
%  - pend_F
%  - pend_H
%

clear variables

% Sampling period
Ts = 0.05;

% Finer discretization
n_step = 5;

% Pendulum parameters
params = struct();
F = 0;
params.g = 10;
params.l = 2;
params.mc = 5;
params.mp = 1;
params.muc = 0.05;
params.mup = 2;
params.dt = Ts / n_step;


%% Compare to DAE simulations

% See file solve_DAE_cartpole.mlx
results_dir = 'results';
filename = "cartpole_benchmark_sim.csv";
DAE_sim_results = readtable( ...
    fullfile(results_dir, filename), ...
    'VariableNamingRule', 'preserve');

%filename = 'pend_sim_benchmark.csv';
%bench_sim_results = readtable(filename,'PreserveVariableNames',true);

% Simulation parameters
t_stop = 20;
N = ceil(t_stop / Ts);  % no. of time steps

% Time series
k = (0:N)';
t = Ts*k;

% Initial system state
var_labels = {'x(t)', 'Dxt(t)', 'theta(t)', 'Dthetat(t)'};
xk = DAE_sim_results{1, var_labels}';

% Specify an input sequence
U = [F*ones(N/2, 1); -F*ones(N/2+1, 1)];
%U(t >= 0.5) = 10;  % Step

% Arrays to save simulation results (u and y)
X = nan(N+1, 4);
Y = nan(N+1, 4);
for i = 1:N+1

    % Save state vector
    X(i, :) = xk';

    % Control input
    uk = U(i, :)';

    % Unmeasured input disturbance
    %pk = bench_sim_results{j, 'p(k)'}';
    pk = 0;

    % Pendulum output measurement
    yk = cartpole_yk(xk, uk+pk, params);

    % Test augmented model functions
    %xak = [xk(1:2); pk];
    %yak = pend_MeasurementFcn(xak, uk, params);

    % Get states at next sample time
    for j = 1:n_step
        xk = cartpole_xkp1(xk, uk+pk, params);
    end
    %xakp1 = pend_StateFcn(xak, uk, params);

    % Save output vector
    Y(i,:) = yk';

end
sim_results = table(k, t, U, X, Y);

% Select main variables
v_sel = [1 3];

% Comparison plot
figure(1); clf
if N >= 50
    marker = '.';
else
    marker = 'o';
end
colorOrder = get(gca, 'ColorOrder');
labels = cell(1, 2 * numel(v_sel));
for i = 1:numel(v_sel)
    plot(t, X(:, v_sel(i)), marker, 'Color', colorOrder(i, :)); 
    hold on
    plot(DAE_sim_results{:, 't'}, DAE_sim_results{:, var_labels(v_sel(i))}, ...
        'LineWidth', 2, 'Color', colorOrder(i, :));
end
ylim([-pi pi])
xlabel("Time ($t$)", 'Interpreter', 'latex')

leg_labels = repmat(var_labels(v_sel), 2, 1) + ["" " DAE"]';
legend(leg_labels, 'Location', 'Best', 'Interpreter', 'latex')
grid on

%TODO: Somethin' ain't right!
% Is DAE wrong or this discrete model?

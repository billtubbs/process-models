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
% To run this test use the following command:
%
% >> runtests test_EKF_observer
%

Ts = 0.1;  % sampling period

% Pendulum parameters
params = struct();
params.F = F;
params.g = g;
params.l = l;
params.mc = mc;
params.mp = mp;
params.muc = muc;
params.mup = mup;
params.dt = Ts;


%% Commpare to DAE simulations

% See file solve_DAE_cartpole.mlx

%filename = 'pend_sim_benchmark.csv';
%bench_sim_results = readtable(filename,'PreserveVariableNames',true);

% Simulation parameters
N = 10;  % no. of time steps

% Time series
k = (0:N)';
t = Ts*k;

% Initial system state
%xk = bench_sim_results{1, {'x1(k+1)', 'x2(k+1)'}}';
xk = [0 0 pi/12 0]';

% Specify an input sequence
U = zeros(N+1, 1);
U(t >= 0.5) = 10;  % Step

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
    xkp1 = cartpole_xkp1(xk, uk+pk, params);
    %xakp1 = pend_StateFcn(xak, uk, params);

    % Save output vector
    Y(i,:) = yk';

    xk = xkp1;

end

sim_results = table(k, t, U, X, Y);
sim_results


% test the following files:
%
% pend_xkp1.m
% pend_yk.m
% pend_StateFcn
% pend_MeasurementFcn.m
% pend_F
% pend_H
%
% To run this test use the following command:
%
% >> runtests test_EKF_observer
%

Ts = 0.1;  % sampling period

% Pendulum parameters
params.K = 1.2;
params.m = 0.3;
params.L = 0.4;
params.g = 9.8;
params.dt = Ts;


%% Example from homework from GEL-7029 course

filename = 'pend_sim_benchmark.csv';
bench_sim_results = readtable(filename,'PreserveVariableNames',true);

% Simulation parameters
N = 150;  % no. of simulation points = 15 sec

% Time series
k = (0:N-1)';
t = Ts*k;

% Initial system state
xk = bench_sim_results{1, {'x1(k+1)', 'x2(k+1)'}}';

% matrix to save the data (u and y)
data = nan(N,11);
for j = 1:N

    % Control input
    uk = bench_sim_results{j, 'u(k)'}';

    % Unmeasured input disturbance
    pk = bench_sim_results{j, 'p(k)'}';

    % Pendulum output measurement
    yk = pend_yk(xk, uk+pk, params);

    % Test augmented model functions
    xak = [xk(1:2); pk];
    yak = pend_MeasurementFcn(xak, uk, params);

    % Get states at next sample time
    xkp1 = pend_xkp1(xk, uk+pk, params);
    xakp1 = pend_StateFcn(xak, uk, params);

    % Record data
    data(j,:) = [k(j) t(j) uk pk xkp1' yk xakp1' yak];

    xk = xkp1;

end

col_names = {'k', 't', 'u(k)', 'p(k)', 'x1(k+1)', 'x2(k+1)', 'y(k)', ...
    'xa1(k+1)', 'xa2(k+1)', 'xa3(k+1)', 'ya(k)'};
sim_results = array2table(data, 'VariableNames', col_names);

% Show selected results
j = find(t == 4.5);
selection = (0:9)' + j;

% Compare states with benchmark data
%sim_results(selection, {'t', 'x1(k+1)', 'x2(k+1)', 'y(k)'})
%bench_sim_results(selection, {'x1(k+1)', 'x2(k+1)', 'y(k)'})
assert(isequal( ...
    round(sim_results{:, {'t', 'x1(k+1)', 'x2(k+1)', 'y(k)'}}, 4), ...
    round(bench_sim_results{:, {'t', 'x1(k+1)', 'x2(k+1)', 'y(k)'}}, 4) ...
))

% Compare augmented model results
%sim_results{selection, {'t', 'xa1(k+1)', 'xa2(k+1)', 'xa3(k+1)', 'ya(k)'}}
%bench_sim_results{selection, {'t', 'x1(k+1)', 'x2(k+1)', 'p(k)', 'y(k)'}}
assert(isequal( ...
   round(sim_results{:, {'t', 'xa1(k+1)', 'xa2(k+1)', 'xa3(k+1)', 'ya(k)'}}, 4), ...
   round(bench_sim_results{:, {'t', 'x1(k+1)', 'x2(k+1)', 'p(k)', 'y(k)'}}, 4) ...
))


%% Test Jacobian calculation functions
%
%   F = df/dxak = pend_F(xak,uk,params)
%   H = dh/dxak = pend_H(xak,uk,params)
%
% where:
%   f = pend_StateFcn
%   h = pend_MeasurementFcn(xak,uk,params)
%

% Number of states of augmented model
na = 3;

% Estimate Jacobian of state transition function
F_est = nan(na, na);
e = 0.0001;
uk = randn();
xak = zeros(na, 1);
for i = 1:na
    for j = 1:na
        dxa = zeros(na, 1);
        dxa(j) = e;
        xak1 = xak - dxa;
        xak2 = xak + dxa;
        f1 = pend_StateFcn(xak1,uk,params);
        f2 = pend_StateFcn(xak2,uk,params);
        F_est(i, j) = (f2(i) - f1(i)) / (2 * dxa(j));
    end
end

%[F_est pend_F(xak,uk,params)]
assert(isequal(round(F_est, 4), round(pend_F(xak,uk,params), 4)))
F_test = [    1.0000    0.1000         0
             -2.4500    0.6000    2.0833
                   0         0    1.0000];
assert(isequal(round(F_est, 4), F_test))

% Estimate Jacobian of measurement function
H_est = nan(1, na);
e = 0.0001;
uk = randn();
xak = zeros(na, 1);
for j = 1:na
    dxa = zeros(na, 1);
    dxa(j) = e;
    xak1 = xak - dxa;
    xak2 = xak + dxa;
    h1 = pend_MeasurementFcn(xak1,uk,params);
    h2 = pend_MeasurementFcn(xak2,uk,params);
    H_est(1, j) = (h2 - h1) / (2 * dxa(j));
end

%[H_est pend_H(xak,uk,params)]
assert(isequal(round(H_est, 4), round(pend_H(xak,uk,params), 4)))
H_test = [1 0 0];
assert(isequal(round(H_est, 4), H_test))

% Tests the following functions:
%
%  - arom3d.m, 
%  - arom3_StateFcn.m
%  - arom3_MeasurementFcn.m
%  - arom3_StateFcnRodin.m
%  - arom3_MeasurementFcnRodin2.m
%  - arom3_StateJacobianFcnRodin.m
%  - arom3_MeasurementJacobianFcnRodin2.m
%
% To run this test use the following command:
%
% >> runtests test_EKF_observer
%

% Initialization for each test
clear all;
seed = 0;
rng(seed)

% Load parameters from file arom3_params.m
arom3_params

% System dimensions
n = size(x0, 1);
nu = size(p0, 1);

y0 = arom3_measurements(x0);
ny = size(y0, 1);

% sample period in minutes
dt = 5/60;


%% Iterative simulation

t_stop = 50;  % length of simulation (hours)
nT = t_stop / dt;  % number of sample periods
t = 0;
x = x0;
p = p0;

X = nan(nT, n);
for k = 1:nT
    [x, y] = arom3d(t, x, p, params, dt);
    X(k,:) = x';
end

% Check steady-state values match those of continuous-time
% model solution (from test_arom3.m)
x0_test = [741.3093;
           466.6142;
           533.3858];
assert(all(abs(X(end,:)' - x0_test) < 0.01))


%% Replicate simulation in Robertson and Lee (1998)

% Note figures 1 and 3 seem to be of the same simulation, not
% figures 1 and 2 as it states in the paper.

dt = 1/60;  % sampling period in minutes
% TODO: Robertson uses different sampling periods for each
% measurement: dt = [1/60; 10/60; 10/60]
t_stop = 15;  % length of simulation (hours)
nT = t_stop / dt;  % duration of simulation

% Measurement noise variance
% (specified as a % of the nominal values)
sigma_v = [0.001; 0.02] .* x0([1 3]);
V = randn(nT, ny) .* sigma_v';

t = 0;
x = x0;
w = zeros(3, 1);

k_ind = (0:nT)';
t_sim = dt*k_ind;

% Try to mimic the disturbance values in Fig. 3
P = zeros(nT+1, 2);
step_periods = {t_sim < 4.6, t_sim >= 4.6, t_sim >= 11.2};
step_values = {[5.4 5.6], [2.7 5.6], [2.7 3]};
for i = 1:numel(step_periods)
    P(step_periods{i}, :) = repmat(step_values{i}, sum(step_periods{i}), 1);
end

X = nan(nT+1, n);
Y_m = nan(nT+1, ny);
x = [740; 479; 532];  % initial condition estimated from fig. 1
for i = 1:nT
    k = k_ind(i);
    p = P(i, :)';
    t = t_sim(k_ind == k);
    X(i, :) = x';
    [xkp1, y] = arom3d(t, x, p, params, dt);
    y_m = y + V(i, :)';
    Y_m(i, :) = y_m';

    % Check observer functions
    xkp1_2 = arom3_StateFcn(x, p, dt, params);
    assert(all(abs(xkp1_2 - xkp1) < 1e-10))
    y_2 = arom3_MeasurementFcn(x, p, dt, params);
    assert(isequal(y_2, y))

    xa = [x; p];
    u = [];
    xakp1 = arom3_StateFcnRodin(xa, u, dt, params);
    assert(all(abs(xakp1 - [xkp1; p]) < 1e-10))

    y_3 = arom3_MeasurementFcnRodin2(xa, u, dt, params);
    assert(isequal(y_3, y))

    x = xkp1;
    
end

% Simulation results
sim_results = table(k_ind, t_sim, P, X, Y_m);

% Values estimated from figure 1 in Robertson et al 1998

T_test = [
      0.0 740.6;
      2.5 737.9;
      5.0 745.0;
      7.5 751.2;
     10.0 750.5;
     12.5 734.4
];

T_sim = nan(size(T_test));
for i = 1:size(T_test, 1)
    t = T_test(i, 1);
    T = T_test(i, 2);
    x_sim = sim_results{t_sim == t, 'X'};
    T_sim(i, :) = [t x_sim(1)];
    %fprintf("%f: %8f %8f\n", t, T, x_sim(1))
end
assert(all(abs(T_sim - T_test) < 1, [1 2]))

% % Plot simulation output
% 
% % Plot evolution of states
% figure(1); clf
% ax1 = subplot(2, 1, 1);
% plot(t_sim, X(:, 1), '-', t_sim, Y_m(:, 1), '.')
% xlabel('t')
% ylabel('Temperature [K]')
% title('Temperature - (x_1)')
% grid on
% 
% ax2 = subplot(2, 1, 2);
% plot(t_sim, X(:, 3), '-', t_sim, Y_m(:, 2), '.')
% xlabel('t')
% ylabel('[gmol/m^3]')
% title('Toluene concentration (x_3)')
% grid on
% 
% linkaxes([ax1 ax2], 'x')
% 
% % Plot disturbance inputs
% figure(2); clf
% ax1 = subplot(2, 1, 1);
% plot(t_sim, P(:, 1))
% xlabel('t')
% ylabel('[h^-1]')
% title('Frequency factor (k_0)')
% grid on
% 
% ax2 = subplot(2, 1, 2);
% plot(t_sim, P(:, 2));
% xlabel('t')
% ylabel('[J/gmol/K]')
% title('Overall heat transfer coefficient (U)')
% grid on
% 
% linkaxes([ax1 ax2], 'x')


%% Test Jacobian calculation functions
%
%   F = df/dxak = arom3_StateJacobianFcnRodin(xak,uk,params)
%   H = dh/dxak = arom3_MeasurementJacobianFcn(xak,uk,params)
%
% where:
%   f = arom3_StateFcnRodin
%   h = arom3_MeasurementJacobianFcnRodin2(xak,uk,params)
%

% Number of states of augmented model
na = 5;

% Estimate Jacobian of state transition function
F_est = nan(na, na);
e = 0.01;
uk = [];
xak = [x0; p0];
for i = 1:na
    for j = 1:na
        dxa = zeros(na, 1);
        dxa(j) = e;
        xak1 = xak - dxa;
        xak2 = xak + dxa;
        f1 = arom3_StateFcnRodin(xak1,uk,dt,params);
        f2 = arom3_StateFcnRodin(xak2,uk,dt,params);
        F_est(i, j) = (f2(i) - f1(i)) / (2 * dxa(j));
    end
end

F_calc = arom3_StateJacobianFcnRodin(xak,uk,dt,params);
% [round(F_est, 4) round(F_calc, 4)]
assert(isequal(round(F_est, 4), round(F_calc, 4)))
F_test = [
    0.7738   -0.0111         0   -1.0237    1.0310;
   -0.1346    0.9819         0   -0.9005         0;
    0.1346    0.0097    0.9917    0.9005         0;
         0         0         0    1.0000         0;
         0         0         0         0    1.0000
];
assert(isequal(round(F_est, 4), F_test))

% Estimate Jacobian of measurement function
H_est = nan(ny, na);
e = 0.0001;
uk = randn();
xak = [x0; p0];
for i = 1:ny
    for j = 1:na
        dxa = zeros(na, 1);
        dxa(j) = e;
        xak1 = xak - dxa;
        xak2 = xak + dxa;
        h1 = arom3_MeasurementFcnRodin2(xak1,uk,dt,params);
        h2 = arom3_MeasurementFcnRodin2(xak2,uk,dt,params);
        H_est(i, j) = (h2(i) - h1(i)) / (2 * dxa(j));
    end
end

H_calc = arom3_MeasurementJacobianFcnRodin2(xak,uk,dt,params);
% [H_est H_calc]

assert(isequal(round(H_est, 4), round(H_calc, 4)))
H_test = [     1     0     0     0     0
               0     0     1     0     0];
assert(isequal(round(H_est, 4), H_test))

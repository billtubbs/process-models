% Test cstr5d.m

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


%% Iterative simulation

dt = 5/60;  % sample period in minutes
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
x = [740; 479; 532];  % estimated from fig. 1 
%x = [737.0; 479.4; 520.6];  % from steady-state simulation
for i = 1:nT
    k = k_ind(i);
    p = P(i, :)';
    t = t_sim(k_ind == k);
    [x, y] = arom3d(t, x, p, params, dt);
    y_m = y + V(i, :)';
    Y_m(i, :) = y_m';
    X(i, :) = x';
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
% Test cstr5d.m

clear all;

% Load parameters from file arom3_params.m
arom3_params

n = size(x0, 1);


%% Iterative simulation

dt = 5/60;  % sample period in minutes
t_stop = 50;  % length of simulation (hours)
nT = t_stop / dt;  % number of sample periods
t = 0;
x = x0;
p = p0;
w = zeros(3, 1);

X = nan(nT, n);
for k = 1:nT
    [x, y] = arom3d(t, x, p, w, params, dt);
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

% Measurement noise
sigma_M = [0.001; 0.02; 0.02] .* x0;  % % of the nominal temperature
V = randn(nT, n) .* sigma_M';

t = 0;
x = x0;
w = zeros(3, 1);

k_ind = (1:nT)';
t_span = dt*k_ind;

% Try to mimic the disturbance values in Fig. 3
P = zeros(nT, 2);
step_periods = {t_span < 4.6, t_span >= 4.6, t_span >= 11.2};
step_values = {[5.4 5.6], [2.7 5.6], [2.7 3]};
for i = 1:numel(step_periods)
    P(step_periods{i}, :) = repmat(step_values{i}, sum(step_periods{i}), 1);
end

X = nan(nT, n);
Y_m = nan(nT, n);
%x = [740; 530; ?];  % estimated from fig. 1 
x = [737.0; 479.4; 520.6];  % from steady-state simulation
for i = 1:nT
    k = k_ind(i);
    p = P(i, :)';
    t = t_span(k_ind == k);
    [x, y] = arom3d(t, x, p, w, params, dt);
    y_m = y + V(i, :)';
    Y_m(i, :) = y_m';
    X(i, :) = x';
end

table(k_ind, t_span, P, X, Y_m)


% Plot evolution of states
figure(1)
ax1 = subplot(3, 1, 1);
plot(t_span, X(:, 1), '-', t_span, Y_m(:, 1), '.')
xlabel('t')
ylabel('Temperature [K]')
title('Temperature - (x_1)')
grid on

ax2 = subplot(3, 1, 2);
plot(t_span, X(:, 2), t_span, Y_m(:, 2), '.')
xlabel('t')
ylabel('[gmol/m^3]')
title('Heptane concentration (x_2)')
grid on

ax3 = subplot(3, 1, 3);
plot(t_span, X(:, 3), t_span, Y_m(:, 3), '.')
xlabel('t')
ylabel('[gmol/m^3]')
title('Toluene concentration (x_3)')
grid on

linkaxes([ax1 ax2 ax3], 'x')


% Plot disturbances
figure(2)
ax1 = subplot(2, 1, 1);
plot(t_span, P(:, 1))
xlabel('t')
ylabel('[h^-1]')
title('Frequency factor (k_0)')
grid on

ax2 = subplot(2, 1, 2);
plot(t_span, P(:, 2));
xlabel('t')
ylabel('[J/gmol/K]')
title('Overall heat transfer coefficient (U)')
grid on

linkaxes([ax1 ax2], 'x')
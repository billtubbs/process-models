% Run a simulation of aromatization process model and save data

clear all; clc;
rng(1,'twister');

%% Prepare directories

data_dir = 'data';
if ~isfolder(data_dir)
    mkdir(data_dir)
end

%% Load system parameters and test system equations

% Load system parameters from file arom3_params.m
arom3_params

n = 3;
ny = 3;
nu = 3;
np = 2;
assert(size(x0, 1) == n)
assert(size(p0, 1) == np)

% Check nominal operating point
t = 0;
x = x0;
w0 = zeros(n,1);  % process disturbances
[dx, y] = arom3(t, x, p0, w0, params);
assert(all(abs(dx) < 0.1))

return


% Check ODE calculation
odefun = @ (t, y) arom3(t, y, p0, w0, params);
t_span = [0 20];
[t, X] = ode45(odefun, t_span, x0);
assert(all(abs(X(end,:)' - x0) < 0.1))


%% Simulate system

% Different simulation types:
sim_names = {'Rob_1', 'Rob_2', 'step_1', 'rand_1'};

% Choose which simulation type to run
sim_case = 3;

fprintf("Simulation case %d : %s\n", sim_case, sim_names{sim_case})

switch sim_case
    case 1  % sim #1 in Robertson

    % Measurement noise variance
    %sigma_v = [50000; 0.05];
    sigma_v = [0; 0]

    % Process noise variance
    w0 = zeros(n,1);  % process disturbances

    % simulation settings
    t_stop = 15;
    nT = t_stop * 60;
    t_sim = linspace(0, t_stop, nT+1)';
    Ts = diff(t_sim(1:2));

    % Set initial input disturbance manually
    p_init = [6.75; 6.6];  % approx values for sim #1 in Robertson

    % find steady-state solution
    odefun = @ (t, y) cstr5(t, y, p_init, w0, params);
    t_span = [0 20];
    [t, X] = ode45(odefun, t_span, x0);
    x_init = X(end,:)'

    % Define input signal p(t) for sim #1 in Robertson
    p_bounds = [4 12;
              2 10];
    p_sim = p_init'.*ones(nT+1, 2);
    p_sim(t_sim >= 2.4, 2) = p_init(2) - 0.5;  % step 1
    p_sim(t_sim >= 10.5, 2) = p_init(2) - 0.3;  % step 2
    p_func = @ (t) interp1(t_sim, p_sim, t, 'previous')';

    case 2  % sim #2 in Robertson
    
    % Measurement noise variance
    %sigma_v = [50000; 0.05];
    sigma_v = [0; 0]

    % Process noise variance
    w0 = zeros(n,1);  % process disturbances

    % simulation settings
    t_stop = 15;
    nT = t_stop * 60;
    t_sim = linspace(0, t_stop, nT+1)';
    Ts = diff(t_sim(1:2));

    % Set initial input disturbance manually
    p_init = [8.4; 6.2];  % approx values for sim #2 in Robertson

    % find steady-state solution
    odefun = @ (t, y) cstr5(t, y, p_init, w0, params);
    t_span = [0 20];
    [t, X] = ode45(odefun, t_span, x0);
    x_init = X(end,:)'

    % Define input signal p(t) for sim #2 in Robertson
    p_bounds = [4 12;
              2 10];
    p_sim = p_init'.*ones(nT+1, 2);
    p_sim(t_sim >= 1, 1) = p_init(1) - 1.5;  % step 1
    p_sim(t_sim >= 3.6, 2) = p_init(2) - 2.5;  % step 2
    p_sim(t_sim >= 4.5, 1) = p_init(1) - 0.5;  % step 3
    p_sim(t_sim >= 4.7, 1) = p_init(1) + 1.6;  % step 4
    p_sim(t_sim >= 6.6, 1) = p_init(1) + 2.4;  % step 5
    p_sim(t_sim >= 6.6, 2) = p_init(2) - 4.1;  % step 6
    p_sim(t_sim >= 10.5, 1) = p_init(1) + 2.6;  % step 7
    p_sim(t_sim >= 14.2, 1) = p_init(1) + 2.7;  % step 8
    p_func = @ (t) interp1(t_sim, p_sim, t, 'previous')';

    case 3  % Simple steps with nominal x0, p0, no noise

    % Process noise variance
    w0 = zeros(n,1);  % process disturbances

    % simulation settings
    t_stop = 15;
    nT = t_stop * 60;
    t_sim = linspace(0, t_stop, nT+1)';
    Ts = diff(t_sim(1:2));

    % Initial parameter values
    p_init = p0;  % p0 = nominal values

    % Measurement noise variance
    sigma_v = [0; 0]
    
    % find steady-state solution
    odefun = @ (t, y) cstr5(t, y, p_init, w0, params);
    t_span = [0 20];
    [t, X] = ode45(odefun, t_span, x0);
    x_init = X(end,:)';

    % Define input signal p(t) for sim #2 in Robertson
    p_bounds = [4 12;
              2 10];
    p_sim = p_init'.*ones(nT+1, 2);
    p_sim(t_sim >= 5, 1) = p_init(1) - 1;  % step 1
    p_sim(t_sim >= 10, 2) = p_init(2) + 1;  % step 2
    p_func = @ (t) interp1(t_sim, p_sim, t, 'previous')';

    case 4  % Pseudo-random (as per Robertson 1995)

    % Measurement noise variance
    %sigma_v = [50000; 0.05];
    sigma_v = [0; 0]

    % Process noise variance
    w0 = zeros(n,1);  % process disturbances

    % Random initial states (see Robertson 1995)
    % x_init = (0.4.*rand(n,1) + 0.8) .* x0;
    % x_init(3) = x0(3);  % except for measured state T

    % simulation settings
    t_stop = 15;
    nT = t_stop * 60;
    t_sim = linspace(0, t_stop, nT+1)';
    Ts = diff(t_sim(1:2));

    % Random initial parameter values (see Robertson 1995)
    p_init = (0.4.*rand(nu,1) + 0.8) .* p0
    
    % find steady-state solution
    odefun = @ (t, y) cstr5(t, y, p_init, w0, params);
    t_span = [0 20];
    [t, X] = ode45(odefun, t_span, x0);
    x_init = X(end,:)';
    
    % Generate random walk disturbances p(t)
    s = RandStream('mt19937ar','Seed',0);  % create a dedicated rng
    epsilon = [0.995; 0.995];
    cov_wp = [1.25; 1.25];
    p_bounds = [4 12;
              2 10];
    p_sim = nan(nT+1, 2);
    p = p_init;
    p_sim(1,:) = p';
    for i=2:nT+1
        gamma = rand(s,2,1) < epsilon;
        wp = sqrt(cov_wp) .* randn(2,1);
        p = p + ((1 - gamma) .* wp);
        p = min(max(p, p_bounds(:,1)), p_bounds(:,2));
        p_sim(i,:) = p';
    end
    p_func = @ (t) interp1(t_sim, p_sim, t, 'previous')';

end


%% Simulate trajectory
odefun = @ (t, x) cstr5(t, x, p_func(t), w0, params);

% Solve initial value problem
[t_sim, X] = ode45(odefun, t_sim, x_init);
nT = size(X,1) - 1;

% Compute outputs - with measurement noise
Y = nan(nT+1, ny);
Ym = nan(size(Y));
for i = 1:nT+1
    Y(i,:) = cstr5_measurements(X(i,:))';
    % Add measurement noise
    vk = randn(2,1).*sqrt(sigma_v);
    Ym(i,:) = Y(i,:) + vk';
end

% Create data table
sim_data = [table(t_sim) array2table(p_sim) array2table(X) array2table(Y) array2table(Ym)];
head(sim_data)

% Save data to file
filename = sprintf('sim_data_%s.csv', sim_names{sim_case});
writetable(sim_data,fullfile(data_dir,filename))


%% Plot simulation output

figure(1); clf
subplot(2,1,1)
stairs(t_sim, p_sim(:,1), 'Linewidth', 2)
ylim(p_bounds(1,:))
xlabel('t (hours)');
ylabel('C_i in (kmol.m^-3)');
grid on
title('Inputs')
subplot(2,1,2)
stairs(t_sim, p_sim(:,2), 'Linewidth', 2)
ylim(p_bounds(2,:))
xlabel('t (hours)');
ylabel('C_m in (kmol.m^-3)');
grid on

figure(2); clf
x_labels = {'C_m','C_I','T','D_0','D_1'};
for i=1:n
    subplot(n,1,i)
    plot(t_sim, X(:,i), 'Linewidth', 2)
    ylabel(x_labels{i});
    grid on
    title(sprintf('x_%d',i))
    if i == n
        xlabel('t (hours)');
    end
end

figure(3); clf
subplot(2,1,1)
plot(t_sim,Y(:,1),'Linewidth',2); hold on
plot(t_sim,Ym(:,1),'.');
xlabel('t (hours)');
ylabel('D_1/D_0');
grid on
title('Outputs')
subplot(2,1,2)
plot(t_sim,Y(:,2),'Linewidth',2); hold on
plot(t_sim,Ym(:,2),'.');
xlabel('t (hours)');
ylabel('T');
grid on
legend('y(k)','y_m(k)')

% Run a simulation of 7-state EKF and save data

clear all; clc;
rng default

%% Load simulation data from file
% Run file test_cstr5_sim.m to generate the data file

data_dir = 'data';
filename = 'sim_data_step_1.csv';
sim_data = readtable(fullfile(data_dir,filename));
fprintf("Simulation data loaded from %s\n", filename);
fprintf("Size: (%d, %d)\n", size(sim_data));
head(sim_data)

%% Load system parameters from file cstr5_params.m
cstr5_params

n = 5;
ny = 2;
nu = 2;
y0 = cstr5_measurements(x0);
assert(size(x0, 1) == n)
assert(size(p0, 1) == nu)
assert(size(y0, 1) == ny)

% Select data to be used for EKF simulation
t_sim = sim_data{:,'t_sim'};
Ts = diff(t_sim(1:2))
% Check sample time is fixed
assert(all((diff(t_sim) - Ts) < 1e-6));
p_sim = sim_data{:,{'p_sim1','p_sim2'}};
nT = size(p_sim,1) - 1;
Ym = sim_data{:,{'Ym1','Ym2'}};
X = sim_data{:,{'X1','X2','X3','X4','X5'}};


%% Simulate EKF filter

% Augmented state vector with input disturbance parameters
xa0 = [x0; p0];  % nominal values
na = size(xa0,1);
xak_est = xa0;  % initial estimate

% Initialize covariance matrix
P = diag([0.1 0.1 1e-3 0.1 0.1 0.1 0.1]);
Q = 0.00625 * eye(na);
R = diag([2e-4 4e-4]);

% % Create MATLAB EKF object
% InitialState = x0;
% EKF = extendedKalmanFilter(@cstr5_dynamics, @cstr5_measurements, InitialState);
% EKF.StateTransitionJacobianFcn = @cstr5_J;
% EKF.MeasurementJacobianFcn = @cstr5_measJac;
% EKF.StateCovariance = P;
% EKF.ProcessNoise = diag(sigma_w);
% EKF.MeasurementNoise = diag(sigma_v);

% Arrays to store results
residuals = nan(nT+1,ny);
%X_cor = nan(nT+1,na);
Xa_est = nan(nT+1,na);
trCovP = nan(nT+1,1);

% Check this matches dt in cstr5_dynamics.m and cstr5_J
assert(abs(Ts - 0.016667) < 1e-4)

for k = 1:nT+1
    
    % Process measurements from simulation above
    yk = Ym(k,:)';
    
    % Actual process input disturbances
    %pk = p_sim(k,:)';
    
    % Process state disturbances
    wk = zeros(n,1);
    
    % Normalize measurements and previous state estimates
    ykn = yk ./ y0;
    xakn_est = xak_est ./ xa0;
    
    % Linearize at current operating point by
    % computing augmented Jacobians F(k) and H(k)
    Fa = cstr7_J(xak_est, wk, params);
    Ha = cstr7_measJac(xak_est);
    
    % Normalize the Jacobians
    Fan = Fa .* xa0' ./ xa0;
    Han = Ha .* xa0' ./ y0;
    
    % Update Kalman Filter gain and covariance matrix
    [K, P] = ekf_update(P, Fan, Han, Q, R);
    
    correction = K * (ykn - Han * xakn_est) .* xa0;
    %correction = K * (yk - Ha * xak_est);
    
    % Compute next state estimates
    xk_est = xak_est(1:n,1);
    pk_est = xak_est(n+1:na,1);
    xak_model = [cstr5_dynamics(xk_est, pk_est, wk, params); pk_est];
    
    xak_est = xak_model + correction;
    Xa_est(k,:) = xak_est';
    
    % Save trace of covariance matrix
    trCovP(k,:) = trace(P);

end

sim_data = [table(t_sim) array2table(X) array2table(p_sim) array2table(Xa_est)];
head(sim_data)
filename = 'sim_data_EKF5.csv';
writetable(sim_data,fullfile(data_dir,filename))


%% Plot EKF simulation results

figure(1); clf
x_labels = {'C_m','C_I','T','D_0','D_1'};
for i=1:n
    subplot(n,1,i)
    plot(t_sim, X(:,i), 'Linewidth', 2); hold on
    plot(t_sim, Xa_est(:,i))
    ylabel(x_labels{i});
    grid on
    title(sprintf('x_%d',i))
    if i == n
        xlabel('t (hours)');
    end
    legend('$\hat{x}(k)$','$x(k)$','Interpreter','Latex')
end

figure(2); clf
subplot(2,1,1)
stairs(t_sim, p_sim(:,1), 'Linewidth', 2); hold on
plot(t_sim, Xa_est(:,n+1))
ylim([6 12])
xlabel('t (hours)');
ylabel('C_i in');
grid on
legend('$Ci_{in}(k)$','$\hat{x}_6(k)$','Interpreter','Latex')
title('Input disturbances')
subplot(2,1,2)
stairs(t_sim, p_sim(:,2), 'Linewidth', 2); hold on
plot(t_sim, Xa_est(:,n+2))
ylim([0 8])
xlabel('t (hours)');
ylabel('C_m in');
grid on
legend('$Cm_{in}(k)$','$\hat{x}_7(k)$','Interpreter','Latex')

figure(3); clf
plot(t_sim,trCovP,'Linewidth',2)
xlabel('t (hours)');
ylabel('tr P(k)');
set(gca, 'YScale', 'log')
grid on
title('Trace of Covariance Matrix')

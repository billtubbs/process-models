% Test Extended Kalman Filter on cstr simulation

clear all;
rng default

% Load parameters from file cstr5_params.m
cstr5_params

n = 5;
ny = 2;
nu = 2;
assert(size(x0, 1) == n)
assert(size(p0, 1) == nu)

% Check nominal operating point
t = 0;
x = x0;
w0 = zeros(n,5);  % process disturbances
[dx, y] = cstr5(t, x, p0, w0, params);
assert(all(abs(dx) < 0.1))

% Measurement noise variance
sigma_v = [20000; 0.2];

% Process noise variance
% TODO: This is not included in the simulation yet
sigma_w = (0.05 * x0).^2;

% Check ODE calculation
odefun = @ (t, y) cstr5(t, y, p0, w0, params);
t_span = [0 20];
[t, X] = ode45(odefun, t_span, x0);
assert(all(abs(X(end,:)' - x0) < 0.1))

% Simulation settings
t_stop = 15;
nT = 15*60;
t_sim = linspace(0, t_stop, nT+1);
Ts = diff(t_sim(1:2));

% Define input signal p(t)
p_sim = p0'.*ones(nT,2);
p_sim(t_sim >= 1, 1) = p0(1) + 1;
p_sim(t_sim >= 5, 2) = p0(2) - 1;
p_sim(t_sim >= 10, 1) = p0(1) - 1;
p_func = @ (t) interp1(t_sim, p_sim, t, 'previous')';

% Simulate trajectory - without disturbances
odefun = @ (t, y) cstr5(t, y, p_func(t), w0, params);
[t, X] = ode45(odefun, t_sim, x0);
nT = size(X,1) - 1;
Y = nan(nT+1, ny);
Ym = nan(size(Y));
for i = 1:nT
    Y(i,:) = cstr5_measurements(X(i,:))';
    % Add measurement noise
    vk = randn(2,1).*sqrt(sigma_v);
    Ym(i,:) = Y(i,:) + vk';
end

%% Plot simulation output

figure(1); clf
subplot(2,1,1)
stairs(t_sim, p_sim(:,1), 'Linewidth', 2)
ylim([6 12])
xlabel('t (hours)');
ylabel('C_i in (kmol.m^-3)');
grid on
title('Inputs')
subplot(2,1,2)
stairs(t_sim, p_sim(:,2), 'Linewidth', 2)
ylim([0 8])
xlabel('t (hours)');
ylabel('C_m in (kmol.m^-3)');
grid on

figure(2); clf
x_labels = {'C_m','C_I','T','D_0','D_1'};
for i=1:n
    subplot(n,1,i)
    plot(t, X(:,i), 'Linewidth', 2)
    ylabel(x_labels{i});
    grid on
    title(sprintf('x_%d',i))
    if i == n
        xlabel('t (hours)');
    end
end


figure(3); clf
subplot(2,1,1)
plot(t,Y(:,1),'Linewidth',2); hold on
plot(t,Ym(:,1),'.');
xlabel('t (hours)');
ylabel('D_1/D_0');
grid on
title('Outputs')
subplot(2,1,2)
plot(t,Y(:,2),'Linewidth',2); hold on
plot(t,Ym(:,2),'.');
xlabel('t (hours)');
ylabel('T');
grid on
legend('y(k)','y_m(k)')


%% State Estimator

InitialState = x0;
EKF = extendedKalmanFilter(@cstr5_dynamics, @cstr5_measurements, InitialState);
EKF.MeasurementJacobianFcn = @cstr5_measJac;
EKF.StateTransitionJacobianFcn = @cstr5_J;
EKF.ProcessNoise = diag(sigma_w);
EKF.MeasurementNoise = diag(sigma_v);

% Simulate EKF
nT = size(X,1) - 1;
assert(size(Y,1) == nT+1);

residBuf = [];
xcorBuf = [];
xpredBuf = [];
trCovP = [];

% Check this matches dt in cstr5_dynamics.m and cstr5_J
assert(abs(Ts - 0.016667) < 1e-4)

for k = 1:nT
    
    % Process measurements from simulation above
    yk = Ym(k,:)';
    
    % Process input disturbances
    pk = p_sim(k,:)';
    
    % Process state disturbances
    wk = zeros(n,1);
    
    % Update Kalman Filter
    [Residual, ResidualCovariance] = residual(EKF, yk);
    [CorrectedState, CorrectedStateCovariance] = correct(EKF, yk);
    [PredictedState, PredictedStateCovariance] = predict(EKF, pk, wk, params);
    
    % Save EKF predictions
    residBuf(k,:) = Residual;
    xcorBuf(k,:) = CorrectedState';
    xpredBuf(k,:) = PredictedState';
    trCovP(k,:) = trace(CorrectedStateCovariance);

end

%% Plot EKF simulation results

figure(4); clf
x_labels = {'C_m','C_I','T','D_0','D_1'};
for i=1:n
    subplot(n,1,i)
    plot(t(2:end), xpredBuf(:,i)); hold on
    plot(t, X(:,i), 'Linewidth', 2)
    ylabel(x_labels{i});
    grid on
    title(sprintf('x_%d',i))
    if i == n
        xlabel('t (hours)');
    end
    legend('$\hat{x}(k)$','$x(k)$','Interpreter','Latex')
end

figure(5); clf
plot(t(2:end),trCovP,'Linewidth',2)
xlabel('t (hours)');
ylabel('tr P(k)');
set(gca, 'YScale', 'log')
grid on
title('Trace of Covariance Matrix')
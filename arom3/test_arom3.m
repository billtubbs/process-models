% Test function arom3.m

clear all;

% Load parameters from file arom3_params.m
arom3_params

n = size(x0, 1);


%% Check nominal operating point

t = 0;
[dx, y] = arom3(t, x0, p0, params);
assert(all(abs(dx ./ x0) < 0.01))

% Actual equilibrium is closer to
% x0 = [741.3090  
%       466.6192;  
%       533.3808];
% See test below


%% Check ODE calculation

odefun = @ (t, x) arom3(t, x, p0, params);
t_span = linspace(0, 100, 6001)';
options = odeset('RelTol', 1e-7);
[~, X] = ode45(odefun, t_span, x0, options);

% % Plot evolution of states
% figure(1)
% ax1 = subplot(3, 1, 1);
% plot(t_span, X(:, 1))
% xlabel('t')
% ylabel('Temperature [K]')
% grid on
% 
% ax2 = subplot(3, 1, 2);
% plot(t_span, X(:, 2));
% xlabel('t')
% ylabel('Concentration [gmol/m^3]')
% grid on
% 
% ax3 = subplot(3, 1, 3);
% plot(t_span, X(:, 3));
% xlabel('t')
% ylabel('Concentration [gmol/m^3]')
% grid on
% 
% linkaxes([ax1 ax2 ax3], 'x')

x0_est = X(end,:)';
assert(all(abs(x0_est - x0) < 5))

% TODO: x2 and x3 are off by 3.16


%% Stability checks

odefun = @ (t, x) arom3(t, x, p0, params);
t_span = linspace(0, 50, 401);
options = odeset('RelTol',1e-6);
dx0 = x0.*0.025;
[t, X] = ode45(odefun, t_span, x0 + dx0, options);
assert(all(abs(X(end,:)' - x0) < 5))
dx0 = -x0.*0.025;
[t, X] = ode45(odefun, t_span, x0 + dx0, options);
assert(all(abs(X(end,:)' - x0) < 5))


%% Part factorial experiment

odefun = @ (t, x) arom3(t, x, p0, params);
t_span = linspace(0, 50, 401);
options = odeset('RelTol',1e-6);
dBB = bbdesign(n,'center',1);
n_exp = length(dBB);
X0 = nan(n_exp,n);
XN = nan(n_exp,n);
settled = nan(n_exp,1);
for i = 1:n_exp
    dx0 = 0.025*x0.*dBB(i,:)';
    X0(i,:) = (x0 + dx0);
    [t, X] = ode45(odefun, t_span, X0(i,:)', options);
    XN(i,:) = X(end,:);
    settled(i) = all(abs(X(end,:)' - x0) < 1e-6);
end

%[table((1:n_exp)', settled) array2table(X0) array2table(XN)]


%% Differences experiment to estimate Jacobian

odefun = @ (t, x) arom3(t, x, p0, params);
J_est = nan(n, n);
e = 0.000001;
for i = 1:n
    for j = 1:n
        dx = zeros(n, 1);
        dx(j) = x0(j) * e;
        x1 = x0 - dx;
        x2 = x0 + dx;
        dxdt1 = odefun(0, x1);
        dxdt2 = odefun(0, x2);
        J_est(i, j) = (dxdt2(i) - dxdt1(i)) / (2 * dx(j));
    end
end

% Test Jacobian
x = x0;
J = arom3_CT_J(x, p0, params);
assert(all(abs(J - J_est) < 1e-10, [1 2]))

J_test = [
   -2.714733   -0.132664         0;
   -1.614948   -0.216695         0;
    1.614948    0.116695   -0.1000];
assert(isequal(round(J, 6), J_test))

% TODO: Need to test arom3_StateJacobianFcn.m

dfdx_est = nan(n, n);
e = 0.0000001;
dt = 1/60;  % Shouldn't make a difference
x = x0;
p = p0;
for i = 1:n
    for j = 1:n
        dx = zeros(n, 1);
        dx(j) = x(j) * e;
        x1 = x - dx;
        x2 = x + dx;
        xkp1_x1 = arom3_StateFcn(x1, p, dt, params);
        xkp1_x2 = arom3_StateFcn(x2, p, dt, params);
        dfdx_est(i, j) = (xkp1_x2(i) - xkp1_x1(i)) ./ (2 * dx(j));
    end
end

% Test Jacobian
dfdx = arom3_StateJacobianFcn(x, p, dt, params);
assert(all(abs(dfdx - dfdx_est) < 5e-2, [1 2]))

%abs(dfdx - dfdx_est)
%TODO: Should this be closer?

dfdx_test = [
    0.9991   -0.0001         0
   -0.0007    0.9999         0
    0.0007    0.0001    1.0000];
assert(isequal(round(dfdx, 4), dfdx_test))


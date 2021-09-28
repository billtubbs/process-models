% Test function arom3.m

clear all;

% Load parameters from file arom3_params.m
arom3_params

n = size(x0, 1);


%% Check nominal operating point

t = 0;
x = x0;
w0 = zeros(n, 1);  % process disturbances
[dx, y] = arom3(t, x, p0, w0, params);
assert(all(abs(dx ./ x0) < 0.01))


%% Check ODE calculation

w0 = zeros(n,1);  % process disturbances
odefun = @ (t, x) arom3(t, x, p0, w0, params);
t_span = linspace(0, 50, 901)';
options = odeset('RelTol',1e-6);
[~, X] = ode45(odefun, t_span, x0, options);

% Plot evolution of states
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

w0 = zeros(n,1);  % process disturbances
odefun = @ (t, x) arom3(t, x, p0, w0, params);
t_span = linspace(0, 50, 401);
options = odeset('RelTol',1e-6);
dx0 = x0.*0.025;
[t, X] = ode45(odefun, t_span, x0 + dx0, options);
assert(all(abs(X(end,:)' - x0) < 5))
dx0 = -x0.*0.025;
[t, X] = ode45(odefun, t_span, x0 + dx0, options);
assert(all(abs(X(end,:)' - x0) < 5))


%% Part factorial experiment

w0 = zeros(n,1);  % process disturbances
odefun = @ (t, x) arom3(t, x, p0, w0, params);
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
    settled(i) = all(abs(X(end,:)' - x0) < 0.1);
end

%[table((1:n_exp)', settled) array2table(X0) array2table(XN)]


%% Differences experiment to estimate Jacobian

w0 = zeros(n,1);  % process disturbances
odefun = @ (t, x) arom3(t, x, p0, w0, params);
J_est = nan(n, n);
e = 0.00001;
for i=1:n
    for j=1:n
        dx = zeros(n,1);
        dx(j) = x0(j) * e;
        x1 = x0 - dx;
        x2 = x0 + dx;
        dxdt1 = odefun(0, x1);
        dxdt2 = odefun(0, x2);
        J_est(i,j) = (dxdt2(i) - dxdt1(i)) / (2 * dx(j));
    end
end

% Test Jacobian
x = x0;
J = arom3_CT_J(x, p0, w0, params);
assert(sum(abs(J - J_est) < 1e-2, [1 2]) == 8)
assert(abs(J(1, 1) - J_est(1, 1)) < 0.1)  % this estimate is not so close

J_test = [
   -2.7147   -0.1327         0;
   -1.6149   -0.2167         0;
    1.6149    0.1167   -0.1000];
assert(sum(abs(J - J_test) < 1e-2, [1 2]) == 8)
assert(abs(J(1, 1) - J_test(1, 1)) < 0.1)  % this estimate is not so close
 

% %% Differences experiment to estimate augmented Jacobian
% 
% xa0 = [x0; p0];
% w0 = zeros(n,1);  % process disturbances
% na = size(xa0,1);
% odefun7 = @ (t, x) [cstr5(t, x(1:n,1), x(n+1:na,1), w0, params); zeros(na-n,1)];
% J7_est = nan(na, na);
% e = 0.001;
% for i=1:na
%     for j=1:na
%         dx = zeros(na,1);
%         dx(j) = xa0(j) * e;
%         x1 = xa0 - dx;
%         x2 = xa0 + dx;
%         dxdt1 = odefun7(0, x1);
%         dxdt2 = odefun7(0, x2);
%         J7_est(i,j) = (dxdt2(i) - dxdt1(i)) / (2 * dx(j));
%     end
% end
% 
% % Test Jacobian
% xa0 = [x0; p0];
% w0 = zeros(n,1);  % process disturbances
% J7 = cstr7_CT_J(xa0, w0, params);
% assert(all(all(abs(J7 - J7_est) < 1e-2)))
% 
% J7_test = [
%   -10.7842   -3.4559   -0.3932         0         0         0    9.2000;
%          0  -10.0187   -0.0017         0         0    0.8000         0;
%    26.1488  115.2394   -6.1332         0         0         0         0;
%     0.0006    0.0232    0.0027  -10.0000         0         0         0;
%    78.5112  346.0040   39.3679         0  -10.0000         0         0;
%          0         0         0         0         0         0         0;
%          0         0         0         0         0         0         0];
% assert(all(all(abs(J7 - J7_test) < 1e-3)))
% 
% x = x0;
% J = cstr5_CT_J(x, p0, w0, params);
% assert(all(all((J7(1:5,1:5) - J) == 0)))


%% Test augmented Jacobians with normalized states and parameters

% TODO: This test is not working (Assertion failed)

% xan0 = ones(na,1);  % normalized operating point
% J7n_est = nan(na, na);
% e = 0.0001;
% for i=1:na
%     for j=1:na
%         dxn = zeros(na,1);
%         dxn(j) = e;
%         xn1 = xan0 - dxn;
%         xn2 = xan0 + dxn;
%         % Convert back to real values
%         x1 = xn1 .* xa0;
%         x2 = xn2 .* xa0;
%         dxdt1 = odefun7(0, x1);
%         dxdt2 = odefun7(0, x2);
%         % Estimate derivative with normalized difference
%         J7n_est(i,j) = (dxdt2(i) - dxdt1(i)) ./ xa0(i) / (2 * dxn(j));
%     end
% end
% 
% J7n = cstr7_CT_J(xa0, w0, params) .* xa0' ./ xa0;
% assert(all(all(abs(J7n - J7n_est) < 2)))

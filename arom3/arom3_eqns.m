% Derive Jacobian symbolically

clear

% Test the function below
assert(strcmp( ...
    replace_array_symbols('y(2) = x1 + x2', 'x', 'dx', 2), ...
    'y(2) = dx(1) + dx(2)'))

% Parameters
syms Ea q V rho Tj Tin A f_star FI Fm Mm mDeltaHp Cp R dt

% Inputs and disturbances
syms Ti Ci_C7H16 Th U k0

% State variables
n = 3;
x = sym('x', [n 1]);
%x = [T;   % reaction temperature, [K]
%     Ch;  % outlet concentration of heptane, [gmol/m^3]
%     Ct]; % outlet concentration of toluene, [gmol/m^3]

% Process inputs (constant here)
u = [Ti;        % inlet temperature (preheated), [K]
     Ci_C7H16;  % inlet concentration of heptane, [gmol/m^3]
     Th];       % heating temperature, [K] (Tj in Robertson et al.)

% ODEs

% dxdt = sym('dx', [n 1]);
% 
% % dT/dt, rate of change of reaction temperature
% % (also x1 in Watanbe and H.)
% dxdt(1) = q / V * (u(1) - x(1)) - DeltaH / (rho * Cp) * rate_const * x(2) ...
%     + U * A / (rho * Cp * V) * (u(3) - x(1));
% 
% % dCh/dt, rate of change of heptane concentratsion
% % (x4 in Watanbe and H.)
% dxdt(2) = q / V * (u(2) - x(2)) - rate_const * x(2);
% 
% % dCt/dt, rate of change of toluene concentration
% % (x2 in Watanbe and H.)
% dxdt(3) = -q / V * x(3) + rate_const * x(2);

% Intermediate calculations

% Heat of reaction [J/gmol] as a function of temperature
% cD = [3.1496e-7 -1.15e-6 -5.536e-2 6.20244e1 2.2026e5];
cD = sym('cD', [1 5]);
T = x(1);
DeltaH = cD * [T^4; T^3; T^2; T; 1];

% Reaction rate constant
rate_const = k0 * exp(-Ea / (R * x(1)));

% State equations
dxdt = sym('dxdt', size(x));

% dT/dt, rate of change of reaction temperature
% (also x1 in Watanbe and H.)
dxdt(1) = q / V * (u(1) - x(1)) - DeltaH / (rho * Cp) * rate_const * x(2) ...
    + U * A / (rho * Cp * V) * (u(3) - x(1));

% dCh/dt, rate of change of heptane concentration
% (x4 in Watanbe and H.)
dxdt(2) = q / V * (u(2) - x(2)) - rate_const * x(2);

% dCt/dt, rate of change of toluene concentration
% (x2 in Watanbe and H.)
dxdt(3) = -q / V * x(3) + rate_const * x(2);

% Compute Jacobian of ODEs
fprintf("\nJacobian:\n")
J = sym('dfdx', [n n]);
for i = 1:n
    for j = 1:n
        J(i, j) = simplify(diff(dxdt(i), x(j)));
        expr = sprintf("%s", J(i, j));
        expr = replace_array_symbols(expr, 'x', 'x', n);
        expr = replace_array_symbols(expr, 'cD', 'cD', 5);
        fprintf("\nJ(%d, %d) = %s;\n", i, j, expr)
    end
end

% Augmented observer model (with unmeasured disturbance
% inputs)
na = n + 2;
xak = sym('xak', [na 1]);

% Continuous-time differential equations for process states
xkp1 = x + dxdt * dt;

% Augmented states are discrete-time integrators
xakp1 = [subs(xkp1, x, xak(1:3)); xak(4:5)];

% Compute Jacobian of state transition function of
% augmented model
fprintf("\nJacobian of augmented StateFcn:\n")
dfdax = sym('dfdax', [na na]);
for i = 1:na
    for j = 1:na
        % Substitute the parameters with the unnormalized state variables
        xakp1B = subs(xakp1(i), [k0, U], [xak(4) .* 1e8, xak(5) .* 1e5]);
        dfdxa(i, j) = simplify(diff(xakp1B, xak(j)));
        expr = sprintf("%s", dfdxa(i, j));
        expr = replace_array_symbols(expr, 'xak', 'xak', na);
        expr = replace_array_symbols(expr, 'cD', 'cD', 5);
        fprintf("\ndfdxa(%d, %d) = %s;\n", i, j, expr)
    end
end

% Function to replace 'x1' etc with 'x(1)' in string
% representation of expression
function newExpr = replace_array_symbols(expr, v, vnew, n)
    newExpr = expr;
    for i = 1:n
        old = sprintf('%s%d', v, i);
        new = sprintf('%s(%d)', vnew, i);
        newExpr = strrep(newExpr, old, new);
    end
end
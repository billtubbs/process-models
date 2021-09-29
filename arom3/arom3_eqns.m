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

% State transition function (for discrete simulation)
xkp1 = x + dt*J;

% Compute Jacobian of state transition function
fprintf("\nJacobian of StateFcn:\n")
dfdx = sym('dfdx', [n n]);
for i = 1:n
    for j = 1:n
        dfdx(i, j) = simplify(diff(xkp1(i), x(j)));
        expr = sprintf("%s", dfdx(i, j));
        expr = replace_array_symbols(expr, 'x', 'xk', n);
        expr = replace_array_symbols(expr, 'cD', 'cD', 5);
        fprintf("\ndfdx(%d, %d) = %s;\n", i, j, expr)
    end
end


function newExpr = replace_array_symbols(expr, v, vnew, n)
    newExpr = expr;
    for i = 1:n
        old = sprintf('%s%d', v, i);
        new = sprintf('%s(%d)', vnew, i);
        newExpr = strrep(newExpr, old, new);
    end
end
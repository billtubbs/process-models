function [x1, y] = cstr5d(t, x0, u, dt)
% [x1, y] = cstr5d(t, x0, u, dt)
% Discrete time approximation of CSTR process model with 
% 5 states
%
% Solves ODE to compute x(k+1).
% 
% See file cstr5.m
%

    % Solve ODE to compute x(k+1)
    odefun = @ (t, x) cstr5(t, x, u);
    t_span = [t t+dt];
    options = odeset('RelTol',1e-6);
    [t, X] = ode45(odefun, t_span, x0, options);
    x1 = X(end,:)';
    [~, y] = cstr5(t, x0, u);

end
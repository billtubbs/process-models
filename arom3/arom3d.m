function [x1, y] = arom3d(t, x0, p, w, params, dt)
% [x1, y] = arom3d(t, x0, p, w, params, dt)
% Discrete time approximation of aromatization process 
% model with 3 states
% 
% Solves ODE to compute x(k+1).
% 
% See file arom3.m
%

    % Solve ODE to compute x(k+1)
    odefun = @ (t, x) arom3(t, x, p, w, params);
    t_span = [t t+dt];
    options = odeset('RelTol', 1e-6);
    [~, X] = ode45(odefun, t_span, x0, options);
    x1 = X(end, :)';
    y = arom3_outputs(x1, p, w);

end
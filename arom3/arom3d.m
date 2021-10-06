function [xkp1, yk] = arom3d(t, xk, pk, params, dt)
% [xkp1, yk] = arom3d(t, xk, pk, params, dt)
% Discrete time approximation of aromatization process 
% model with 3 states.
% 
% Solves ODE to compute x(k+1).
% 
% See file arom3.m
%

    % Solve ODE to compute x(k+1)
    odefun = @ (t, x) arom3(t, xk, pk, params);
    t_span = [t t+dt];
    options = odeset('RelTol', 1e-6);
    [~, X] = ode45(odefun, t_span, xk, options);
    xkp1 = X(end, :)';
    yk = arom3_outputs(xk, pk);

end
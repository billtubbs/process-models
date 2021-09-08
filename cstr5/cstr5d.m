function [x1, y] = cstr5d(t, x0, u, dt)
    % Discrete time approximation of CSTR model based on :
    % - Daoutidis, Soroush, and Kravaris, 1990 and
    %   presented by Robertson, Kesavan and Lee, 1995
    % 
    % See file cstr5.m
    %
    % Solve ODE to compute x(k+1)
    odefun = @ (t, x) cstr5(t, x, u);
    t_span = [t t+dt];
    [t X] = ode45(odefun, t_span, x0);
    x1 = X(end,:)';
    [~, y] = cstr5(t, x0, u);

end
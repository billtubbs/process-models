%% Jacobian of measurement equations for nonlinear model of CSTR reactor

function H = cstr5_measJac(x)
    % Jacobian of measurement function for the CSTR model
    %
    % See :
    % - Daoutidis, Soroush, and Kravaris, 1990, and
    % - Robertson, Kesavan and Lee, 1995
    %
    % Measurements
    % y(1) = [x(5) / x(4); ... % D1/D0 the number average molecular weight
    % y(2) = x(3)];        ... % Reactor temperature.
    
    % Jacobian of meaurement equations
    H = [0   0   0  -x(5)*x(4)^-2  1/x(4);
         0   0   1              0       0];

end
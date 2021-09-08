%% Jacobian of measurement equations for nonlinear model of CSTR reactor

function Ha = cstr7_measJac(x)
    % Jacobian matrix of measurement function for
    % augmented discrete time CSTR model
    %
    % See :
    % - Daoutidis, Soroush, and Kravaris, 1990, and
    % - Robertson, Kesavan and Lee, 1995
    %
    % Measurements
    % y(1) = [x(5) / x(4); ... % D1/D0 the number average molecular weight
    % y(2) = x(3)];        ... % Reactor temperature.
    
    % Jacobian of meaurement equations
    Ha = [0   0   0  -x(5)*x(4)^-2  1/x(4)   0   0;
          0   0   1              0       0   0   0];

end
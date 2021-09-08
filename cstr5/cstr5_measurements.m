%% Measurements equations for nonlinear model of CSTR reactor

function y = cstr5_measurements(x)
    % CSTR model
    %
    % See :
    % - Daoutidis, Soroush, and Kravaris, 1990, and
    % - Robertson, Kesavan and Lee, 1995
    
    % Output equations
    y = [x(5)/x(4);   ... % D1/D0 the number average molecular weight
              x(3)];  ... % Reactor temperature.

end
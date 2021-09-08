%% Output equation for nonlinear model of CSTR reactor

function y = cstr5_outputs(x, p, w)
    % CSTR model
    %
    % See :
    % - Daoutidis, Soroush, and Kravaris, 1990, and
    % - Robertson, Kesavan and Lee, 1995
    
    % Output equations (outputs = measurements)
    y = cstr5_measurements(x);

end
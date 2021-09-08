function F = cstr5_J(x, u, w, params)
    % Jacobian matrix of discrete time CSTR model based on :
    % - Daoutidis, Soroush, and Kravaris, 1990 and
    %   presented by Robertson, Kesavan and Lee, 1995

    % Sample time (hours)
    dt = 0.016667; 
    
    % Discrete-time Jacobian of state transition function
    F = eye(5) + dt * cstr5_CT_J(x, u, w, params);
    
end
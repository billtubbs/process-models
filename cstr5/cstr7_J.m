function Fa = cstr7_J(x, w, params)
    % Jacobian matrix of augmented discrete time CSTR model 
    % for EKF state estimation based on :
    % - Daoutidis, Soroush, and Kravaris, 1990 and
    %   presented by Robertson, Kesavan and Lee, 1995

    % Sample time (hours)
    dt = 0.016667; 
    
    % Discrete-time Jacobian of state transition function
    Fa = eye(7) + dt * cstr7_CT_J(x, w, params);
    
end
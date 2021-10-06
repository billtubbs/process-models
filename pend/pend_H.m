function H = pend_H(xk,uk,params)
% H = pend_H(xk,uk,params)
% Jacobian of the output measurement function
% for linearizing the pendulum system.
    % For testing ekf_update.m
    assert(isequal(fieldnames(params), {'K', 'm', 'L', 'g', 'dt'}'))
    H = [1 0 1];
end
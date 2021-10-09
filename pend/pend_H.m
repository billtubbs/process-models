function H = pend_H(xak,uk,params)
% H = pend_H(xk,uk,params)
% Jacobian of the output measurement function
% for linearizing the augmented pendulum system.
    % For testing ekf_update.m
    assert(isequal(fieldnames(params), {'K', 'm', 'L', 'g', 'dt'}'))

    H = [1 0 0];

end
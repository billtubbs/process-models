function yk = pend_yk(xk,uk,params)
% yk = pend_yk(xk,uk,params)
% Measurement equation for the non-linear pendulum 
% system.
%
% Arguments:
% xk : state vector [x1; x2] where
%     xk(1) : angle
%     xk(2) : angular velocity
% uk : torque (not used by this function).
% params : struct (not used by this function).
%
    yk = xk(1);
end
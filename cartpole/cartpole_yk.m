function yk = cartpole_yk(xk,uk,params)
% yk = cartpole_yk(xk,uk,params)
% Measurement equation for the non-linear pendulum 
% system.
%
% Arguments:
%   xk : (4, 1) double vector
%       System states at time k, defined:
%           xk(1) : position of cart, x(t), from left to right
%           xk(2) : velocity of cart, Dx(t), from left to right
%           xk(3) : angle of pole, theta(t), clockwise from vertical
%           xk(4) : angular velocity of pole, Dtheta(t), clockwise.
%   uk : (1, 1) double
%       Input variable : force acting on cart (+ve to right).
%   params : struct
%       This is not used by this function.
%

    % For unit tests
    assert(isequal(size(xk), [4 1]))
    assert(isscalar(uk))

    % y(k)
    yk = xk;

end
function xkp1 = pend_xkp1(xk,uk,params)
% xkp1 = pend_xkp1(xk,uk,params)
% State transition function for the non-linear
% pendulum system.
%
% Arguments:
%   xk : state vector [x1; x2] where
%     xk(1) : angle
%     xk(2) : angular velocity.
%   uk : Torque input.
%   params : struct containing the following 
%     parameter values:
%       params.K : coefficient of friction
%       params.m : mass
%       params.L : length
%       params.g : gravitational acceleration
%       params.dt : sampling interval
%
    assert(isequal(size(xk), [2 1]))
    assert(isscalar(uk))
    K = params.K;
	m = params.m;
	L = params.L;
	g = params.g;
	dt = params.dt;

    % x(k+1)
    xkp1 = xk + dt*[xk(2); 
                    -(g/L*sin(xk(1)) + K/m*xk(2))] ...
           + [0; dt/(m*L^2)] * uk;

end
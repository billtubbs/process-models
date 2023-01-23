function xkp1 = cartpole_xkp1(xk,uk,params)
% xkp1 = cartpole_xkp1(xk,uk,params)
% State transition function for simulating a 
% discrete-time model of the non-linear cart-pole 
% system using the DAEs defined by Florian (2007).
%
% See script 'cartpole_DAEs.m'.
%
% Arguments:
%   xk : (5, 1) double vector
%       System states at time k, defined:
%           xk(1) : position of cart, x(t), from left to right
%           xk(2) : velocity of cart, Dx(t), from left to right
%           xk(3) : angle of pole, theta(t), clockwise from vertical
%           xk(4) : angular velocity of pole, Dtheta(t), clockwise.
%   uk : (1, 1) double
%       Input variable : force acting on cart (+ve to right).
%   params : struct
%       Parameter values defined:
%           params.g : acceleration due to gravity (+ve down)
%           params.l : pendulum length
%           params.mc : cart mass
%           params.mp : pole mass
%           params.muc : coefficient of friction for cart and track 
%           params.mup : coefficient of friction for pole joint
%           params.dt : sampling interval.
%

    assert(isequal(size(xk),[4 1]))
    assert(isscalar(uk))

    % Sampling interval
    dt = params.dt;

    % Set F parameter based on input variableParameter values
    params.F = uk;

    % DAE system variables defined as follows:
    %   Y : (5, 1) double vector
    %     Y(1) : horizontal position of cart, x(t), from left to right
    %     Y(2) : angle of pole, theta(t), clockwise from vertical
    %     Y(3) : downward force of cart on track, Nc(t)
    %     Y(4) : velocity of cart, Dxt(t), from left to right
    %     Y(5) : angular velocity of pole, Dthetat(t), clockwise.

    % Initial condition
    y0 = [ ...
        xk(1);  % x(t)
        xk(3);  % theta(t)
        nan;    % Nc(t) is not used by cartpole_DAEs
        xk(2);  % Dx(t)
        xk(4)   % Dtheta(t)
    ];
    % TODO: Is this correct?
    yp0 = zeros(5, 1);

    % Define function handle for DAEs with parameter 
    % values provided
    F_DAE = @(t, Y, YP) cartpole_DAEs(t, Y, YP, params);

    % Determine complete initial condition numerically
    FIXED_Y0 = [1 1 0 1 1]';  % specify which are fixed
    FIXED_YP0 = [0 0 0 1 1]';
    [y0, yp0] = decic(F_DAE, 0, y0, FIXED_Y0, yp0, FIXED_YP0, opt);

    % Solve DAEs with variable step, variable order solver
    opt = odeset('RelTol', 10.0^(-7), 'AbsTol', 10.0^(-7));
    [tSol, ySol] = ode15i(F_DAE, [0 dt], y0, yp0, opt);

    % x(k+1)
    xkp1 = ySol;

end
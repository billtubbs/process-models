function xakp1 = pend_StateFcn(xak,uk,params)
% xakp1 = pend_StateFcn(xak,uk,params)
% State transition function for the non-linear pendulum 
% system augmented with an output disturbance.
%
% Arguments:
% xak : state vector [x1; x2] where
%     xak(1) : angle
%     xak(2) : angular velocity
%     xak(3) : input disturbance model state.
% uk : Torque input.
% params : struct containing the parameter values
%     listed below.
%
% params.K : coefficient of friction
% params.m : mass
% params.L : length
% params.g : gravitational acceleration
% params.T : sample period
%
    % For unit tests
    assert(isequal(fieldnames(params), {'K', 'm', 'L', 'g', 'dt'}'))
    assert(isequal(size(xak),[3 1]))
    assert(isscalar(uk))

    xakp1 = [pend_xkp1(xak(1:2),uk+xak(3),params); xak(3)];

end
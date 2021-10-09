function y = pend_MeasurementFcn(xak,uk,params)
% y = pend_MeasurementFcn(xak,uk,params)
% Measurement equation for the non-linear pendulum 
% system augmented with an output disturbance.
%
% Arguments:
% xa : state vector [x1; x2] where
%     x1 : angle
%     x2 : angular velocity
%     xi : integrator state.
% u : torque (not used by this function).
% params : struct (not used by this function).
%
    % For unit tests
    % For unit tests
    assert(isequal(fieldnames(params), {'K', 'm', 'L', 'g', 'dt'}'))
    assert(isequal(size(xak),[3 1]))
    assert(isscalar(uk))

    y = xak(1);

end
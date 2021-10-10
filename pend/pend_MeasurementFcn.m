function y = pend_MeasurementFcn(xak,uk,params)
% y = pend_MeasurementFcn(xak,uk,params)
% Measurement equation for the non-linear pendulum 
% system augmented with an output disturbance.
%
% Arguments:
% xak : augmented state vector [xak1; xak2; xak3] where
%     xak(1) : angle
%     xak(2) : angular velocity
%     xak(3) : input disturbance model state.
% uk : torque input (not used by this function).
% params : struct (not used by this function).
%
    % For unit tests
    % For unit tests
    assert(isequal(fieldnames(params), {'K', 'm', 'L', 'g', 'dt'}'))
    assert(isequal(size(xak),[3 1]))
    assert(isscalar(uk))

    y = xak(1);

end
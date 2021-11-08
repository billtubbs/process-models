function H = arom3_MeasurementJacobianFcnRodin1(xa, uk, dt, params)
% H = arom3_MeasurementJacobianFcnRodin1(xa, uk, dt, params)
% Jacobian matrix of measurement function of aromatization
% process with 3 states and 1 measured output (T).
%
% State variables
% xa(1) : T, reaction temperature, [K]
% xa(2) : Ch, outlet concentration of heptane, [gmol/m^3]
% xa(3) : Ct, outlet concentration of toluene, [gmol/m^3]
%
% Random shock signals to generate RODD disturbances
% xa(4) : Wp(1), RODD step disturbance applied to k0, the 
%     frequency factor, or pre-exponential factor [1/h].
% xa(5) : Wp(2), RODD step disturbance applied to U, the
%     overall heat transfer coefficient [J/(gmol.K)]
% 
% Output measurements
% y(1) : x(1)
% y(2) : x(3)

    % Jacobian of measurement equations
    H = [1   0   0   0   0];

end
function Ha = arom3_MeasurementJacobianFcnRodin1(xa)
% Ha = arom3_MeasurementJacobianFcnRodin1(xa)
% Jacobian matrix of measurement function of aromatization
% process with 3 states and 1 measured output (T).
%
% State variables
% xa(1) : T, reaction temperature, [K]
% xa(2) : Ch, outlet concentration of heptane, [gmol/m^3]
% xa(3) : Ct, outlet concentration of toluene, [gmol/m^3]
%
% Process disturbances (time-varying parameters)
% xa(4) : k0, frequency factor, or pre-exponential factor [1/h]
% xa(5) : U, overall heat transfer coefficient [J/(gmol.K)]
% 
% Output measurements
% y(1) : x(1)
% y(2) : x(3)

    % Jacobian of meaurement equations
    Ha = [1   0   0   0   0];

end
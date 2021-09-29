function Ha = arom3_MeasurementJacobianFcn(x)
% Ha = arom3_MeasurementJacobianFcn(x)
% Jacobian matrix of measurement function of aromatization
% process with 3 states
%
% State variables
% x(1) : T, reaction temperature, [K]
% x(2) : Ch, outlet concentration of heptane, [gmol/m^3]
% x(3) : Ct, outlet concentration of toluene, [gmol/m^3]
%
% Process inputs (disturbances)
% p(1) : k0, frequency factor, or pre-exponential factor [1/h]
% p(2) : U, overall heat transfer coefficient [J/(gmol.K)]
% 
% Output measurements
% y(1) : x(1)
% y(2) : x(3)

    % Jacobian of meaurement equations
    Ha = [1   0   0;
          0   0   1];

end
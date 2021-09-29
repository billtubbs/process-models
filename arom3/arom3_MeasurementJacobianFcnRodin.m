function Ha = arom3_measJacRodin(xa)
% Ha = arom3_measJacRodin(xa)
% Jacobian matrix of measurement function of aromatization
% process augmented with two input disturbances.
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
    Ha = [1   0   0   0   0;
          0   0   1   0   0];

end
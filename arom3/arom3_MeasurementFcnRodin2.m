function y = arom3_MeasurementFcnRodin2(xa)
% y = arom3_MeasurementFcnRodin2(xa)
% Measurement equations for aromatization process augmented
% with two input disturbances.
%
% Description:
% - Heptane to toluene aromatization process model described
%   by Watanabe and Himmelblau (1983).
%
% References:
% - Watanabe, K. and D. M. Himmelblau (1983). Fault diagnosis
%   in nonlinear chemical processes. AIChE J. 29, 243-261.
%     Part I. Theory
%     Part II. Application to a Chemical Reactor
% - Robertson, D. G., and Lee, J. H. (1998). A Method for the
%   Estimation of Infrequent Abrupt Changes in Nonlinear
%   Systems.
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
%
% Nominal values for the states and parameters
% x0 = [742.0;  % reaction temperature, [K]
%       463.0;  % outlet concentration of heptane, [gmol/m^3]
%       537.0];  % outlet concentration of toluene, [gmol/m^3]
% p0 = [5e8 / 1e8;  % (normalized)
%       6e5 / 1e5];  % (normalized)

    % Output equations
    y = [xa(1);  % T, reaction temperature, [K]
         xa(3)];  % Ct, outlet concentration of toluene, [gmol/m^3]

end
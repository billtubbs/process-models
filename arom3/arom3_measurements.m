function y = arom3_measurements(x)
% y = arom3_measurements(x)
% Measurement equations for aromatization process with 3 states
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
%
% Nominal values for the states and parameters
% x0 = [742.0;  % reaction temperature, [K]
%       463.0;  % outlet concentration of heptane, [gmol/m^3]
%       537.0];  % outlet concentration of toluene, [gmol/m^3]
% p0 = [5e8 / 1e8;  % (normalized)
%       6e5 / 1e5];  % (normalized)

    % Output equations
    y = [x(1);  % T, reaction temperature, [K]
         x(3)];  % Ct, outlet concentration of toluene, [gmol/m^3]

end
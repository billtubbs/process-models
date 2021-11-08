function y = arom3_MeasurementFcnRodin1(xa, uk, dt, params)
% y = arom3_MeasurementFcnRodin1(xa, uk, dt, params)
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
% Random shock signals to generate RODD disturbances
% xa(4) : Wp(1), RODD step disturbance applied to k0, the 
%     frequency factor, or pre-exponential factor [1/h].
% xa(5) : Wp(2), RODD step disturbance applied to U, the
%     overall heat transfer coefficient [J/(gmol.K)]
% 
% Output measurements
% y(1) : x(1)
%
% Nominal values for the states and parameters
% x0 = [742.0;  % reaction temperature, [K]
%       463.0;  % outlet concentration of heptane, [gmol/m^3]
%       537.0];  % outlet concentration of toluene, [gmol/m^3]
% p0 = [5e8 / 1e8;  % (normalized)
%       6e5 / 1e5];  % (normalized)

    % Output equations
    y = xa(1);  % T, reaction temperature, [K]

end
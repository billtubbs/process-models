function xakp1 = arom3_StateFcnRodin(xak, uk, dt, params)
% xakp1 = arom3_StateFcnRodin(xak, uk, dt, params)
% State transition function for the discrete-time model of 
% the aromatization process augmented with two input
% disturbances.
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
% y(2) : x(3)
%
% Nominal values for the states and parameters
% xa0 = [742.0;  % reaction temperature, [K]
%        463.0;  % outlet concentration of heptane, [gmol/m^3]
%        537.0;  % outlet concentration of toluene, [gmol/m^3]
%        5e8 / 1e8;  % p(1) (normalized)
%        6e5 / 1e5];  % p(2) (normalized)

    % Calculate next states of process model;
    xk = xak(1:3);
    pk = p0 + xak(4:5);
    xkp1 = arom3_StateFcn(xk, pk, dt, params);

    % Augmented states are discrete-time integrators
    xakp1 = [xkp1; xak(4:5) + uk];

end
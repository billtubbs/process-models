% Constants and parameters for the heptane to toluene
% aromatization process model presented by Watanabe and 
% Himmelblau (1983) and used by Robertson and Lee (1998)
%
% References:
% - Watanabe, K. and D. M. Himmelblau (1983). Fault diagnosis
%   in nonlinear chemical processes. AIChE J. 29, 243-261.
%     Part I. Theory
%     Part II. Application to a Chemical Reactor
% - Robertson, D. G., and Lee, J. H. (1998). A Method for the
%   Estimation of Infrequent Abrupt Changes in Nonlinear
%   Systems.

% Chemical symbols:
% Heptane, C7H16
% Toluene, C7H8

% Constant parameter values from Watanabe and
% Himmelblau (1983):
params.k0 = 5.01e8;  % reaction rate coefficient (1/h)
params.Ea = 1.369e5;  % activation energy (J/gmol)
params.R = 8.319;  % gas constant (J/(gmol.K))

% The Heat of reaction DeltaH and specific heat Cp
% are given as functions of temperature but only DeltaH
% is considered as a function of temperature because
% the change of Cp with temperature is slight in the
% range of interest here:

% Heat of reaction [J/gmol] as a function of temperature
params.DeltaH = @(T) 2.2026e5 + 6.20244e1*T - 5.536e-2*T^2 - 1.15e-6*T^3 + 3.1496e-7*T^4;

% Specific heat (molar)
params.Cp = 490.7;  % [J/(gmol.K)]

% Other process parameters are assumed to be constant
params.rho = 593.0;  % density, [gmol/m3]
params.h = 6.069e5;  % overall heat transfer coefficient, [J/(gmol.K)]
params.A = 10.0;  % area of heat exchange, [m^2]
params.V = 30.0;  % effective reactor volume, [m^3]
params.q = 3.0;  % inlet and outlet volumetric flow rate, [m^3/h]

% Note: In Watanabe and H., inlet temperature, heating temperature,
% inlet heptane concentration, and outlet H2 concentration, are 
% process variables but in Robertson et al. they are constants:
params.Ti = 600.0;  % inlet temperature (preheated), [K]
params.Th = 850.0;  % heating temperature, [K]
params.Ci_C7H16 = 1000.0;  % inlet concentration of heptane, [gmol/m^3]
params.Co_H2 = 2149.0;  % outlet concentration of H2, [gmol/m^3]

% Process variables at steady-state
params.T0 = 742.0;  % reaction temperature, [K]
params.Co_C7H16 = 463.0;   % outlet concentration of heptane, [gmol/m^3]
params.Co_C7H8 = 537.0;   % outlet concentration of toluene, [gmol/m^3]

% Nominal values of parameters
params.k0 = 5e8;  % frequency factor, or pre-exponential factor [h^-1]
                  % (Watanabe and H. use 5.01e8)
params.U = 6e5;  % overall heat transfer coefficient [J/(gmol.K)]

% Nominal values of the states
x0 = [params.T0;   % reaction temperature, [K]
      params.Co_C7H16;   % outlet concentration of heptane, [gmol/m^3]
      params.Co_C7H8];   % outlet concentration of toluene, [gmol/m^3]

% TODO: Currently, actual equilibrium point of the model is
% closer to the following check equations to see if there 
% is an error
% x0 = [741.3094
%       466.6194
%       533.3806];

% Normalized nominal input values (disturbances)
p0 = [params.k0 / 1e8;
      params.U / 1e5];


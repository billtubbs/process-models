function dx = arom3_dynamics_CT(x, p, w, params)
% dx = arom3_dynamics_CT(x, p, w, params)
% Continous-time model of aromatization process with 3 states
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
% Process disturbances
% p(1) : k0, frequency factor, or pre-exponential factor [1/h]
% p(2) : U, overall heat transfer coefficient [J/(gmol.K)]
% 
% Nominal values for the states and parameters
% x0 = [742.0;  % reaction temperature, [K]
%       463.0;  % outlet concentration of heptane, [gmol/m^3]
%       537.0];  % outlet concentration of toluene, [gmol/m^3]
% p0 = [5e8 / 1e8;  % (normalized)
%       6e5 / 1e5];  % (normalized)
%
% See documentation on defining grey-box model files:
% https://www.mathworks.com/help/ident/ug/creating-idnlgrey-model-files.html

    % Constant parameter values:
    k0 = params.k0;  % frequency factor (1/h)
    Ea = params.Ea;  % activation energy (J/gmol)
    R = params.R;  % gas constant (J/(gmol.K))

    % Heat of reaction [J/gmol] as a function of temperature
    DeltaH = params.DeltaH(x(1));

    % Specific heat (molar)
    Cp = params.Cp;  % [J/(gmol.K)]

    % Other process parameters are assumed to be constant
    rho = params.rho;  % density, [gmol/m3]
    h = params.h;  % overall heat transfer coefficient, [J/(gmol.K)]
    A = params.A;  % area of heat exchange, [m^2]
    V = params.V;  % effective reactor volume, [m^3]
    q = params.q;  % inlet and outlet volumetric flow rate, [m^3/h]

    % Nominal values of the states and parameters
    x0 = [params.T0;   % reaction temperature, [K]
          params.Co_C7H16;   % outlet concentration of heptane, [gmol/m^3]
          params.Co_C7H8];   % outlet concentration of toluene, [gmol/m^3]

    % Parameter values from normalized values in p(k)
    k0 = p(1) .* 1e8;  % frequency or pre-exponential factor [h^-1]
    U = p(2) .* 1e5;  % overall heat transfer coefficient [J/(gmol.K)]

    % Process inputs (constant here)
    u = [params.Ti;  % inlet temperature (preheated), [K]
         params.Ci_C7H16;  % inlet concentration of heptane, [gmol/m^3]
         params.Th];  % heating temperature, [K] (Tj in Robertson et al.)
    
    % Add process disturbances to states
    u = u + w;

    % Reaction rate constant
    rate_const = k0 * exp(-Ea / (R * x(1)));

    % State equations
    dx = nan(3, 1);

    % dT/dt, rate of change of reaction temperature
    % (also x1 in Watanbe and H.)
    dx(1) = q / V * (u(1) - x(1)) - DeltaH / (rho * Cp) * rate_const * x(2) ...
        + U * A / (rho * Cp * V) * (u(3) - x(1));

    % dCh/dt, rate of change of heptane concentration
    % (x4 in Watanbe and H.)
    dx(2) = q / V * (u(2) - x(2)) - rate_const * x(2);

    % dCt/dt, rate of change of toluene concentration
    % (x2 in Watanbe and H.)
    dx(3) = -q / V * x(3) + rate_const * x(2);

end
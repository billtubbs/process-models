%% Output equation for nonlinear model of aromatization process

function y = arom3_outputs(x, p, w)
    % Aromatization process model
    %
    % Description:
    % - Heptane to toluene aromatization process described by
    %   Watanabe and Himmelblau (1983).
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
    
    % Output equations (outputs = measurements)
    y = arom3_measurements(x);

end
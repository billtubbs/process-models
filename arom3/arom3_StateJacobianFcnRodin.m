function dfdxa = arom3_StateJacobianFcnRodin(xak, dt, params)
% dfdxa = arom3_StateJacobianFcnRodin(xak, dt, params)
% Jacobian of the state transition function for the discrete-
% time model of the aromatization process augmented with two 
% input disturbances.
%
% Description:
% - Heptane to toluene aromatization process model described
%   by Watanabe and Himmelblau (1983).
%
% References:
% - Watanabe, K. and D. M. Himmelblau (1983). Fault diagnosis
%   in nonlinear chemical processes. AIChE J. 29, 243-261.
%     Part I. xak(1)heory
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
% y(1) : xak(1)
% y(2) : xak(3)
%
% Nominal values for the states and parameters
% x0 = [742.0;  % reaction temperature, [K]
%       463.0;  % outlet concentration of heptane, [gmol/m^3]
%       537.0];  % outlet concentration of toluene, [gmol/m^3]
% p0 = [5e8 / 1e8;  % (normalized)
%       6e5 / 1e5];  % (normalized)

    % Parameter values from normalized values in p(k)
    k0 = xak(4) .* 1e8;  % frequency or pre-exponential factor [h^-1]
    U = xak(5) .* 1e5;  % overall heat transfer coefficient [J/(gmol.K)]

    % Constant parameter values:
    Ea = params.Ea;  % activation energy (J/gmol)
    R = params.R;  % gas constant (J/(gmol.K))

    % Specific heat (molar)
    Cp = params.Cp;  % [J/(gmol.K)]

    % Other process parameters are assumed to be constant
    rho = params.rho;  % density, [gmol/m3]
    h = params.h;  % overall heat transfer coefficient, [J/(gmol.K)]
    A = params.A;  % area of heat exchange, [m^2]
    V = params.V;  % effective reactor volume, [m^3]
    q = params.q;  % inlet and outlet volumetric flow rate, [m^3/h]

    % Heat of reaction [J/gmol] as a function of temperature
    cD = params.cD;
    DeltaH = cD * [xak(1)^4; xak(1)^3; xak(1)^2; xak(1); 1];

    % Partial derivative
    dDeltaHdx1 = (cD(4) + 2*cD(3)*xak(1) + 4*cD(1)*xak(1)^3 + 3*cD(2)*xak(1)^2);

    % Reaction rate constant
    rate_const = k0 * exp(-Ea / (R * xak(1)));

    % Jacobian matrix
    dfdxa = zeros(5, 5);

    dfdxa(1, 1) = 1 - dt*( ...
            (xak(2)*rate_const*(2*cD(3) + 6*cD(2)*xak(1) + 12*cD(1)*xak(1)^2)) ...
            / (Cp*rho) + (Ea^2*xak(2)*rate_const*DeltaH) / (Cp*R^2*rho*xak(1)^4) ...
            + (2*Ea*xak(2)*rate_const*dDeltaHdx1)/(Cp*R*rho*xak(1)^2) ...
            - (2*Ea*xak(2)*rate_const*DeltaH)/(Cp*R*rho*xak(1)^3) ...
        );

    dfdxa(1, 2) = -dt*( ...
            (rate_const*dDeltaHdx1)/(Cp*rho) ...
            + (Ea*rate_const*DeltaH)/(Cp*R*rho*xak(1)^2) ...
        );

    dfdxa(2, 1) = -(Ea*dt*xak(2)*rate_const*(Ea - 2*R*xak(1)))/(R^2*xak(1)^4);

    dfdxa(2, 2) = 1 - (Ea*dt*rate_const)/(R*xak(1)^2);

    dfdxa(3, 1) = (Ea*dt*xak(2)*rate_const*(Ea - 2*R*xak(1)))/(R^2*xak(1)^4);

    dfdxa(3, 2) = (Ea*dt*rate_const)/(R*xak(1)^2);

    % TODO: more derivative terms to add
    
    dfdxa(3, 3) = 1;

    dfdxa(4, 4) = 1;

    dfdxa(5, 5) = 1;

end
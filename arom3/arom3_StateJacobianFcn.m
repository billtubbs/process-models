function dfdx = arom3_StateJacobianFcn(xk, pk, dt, params)
% dfdx = arom3_StateJacobianFcn(xk, pk, dt, params)
% Jacobian of the state transition function for the  
% discrete-time model of the aromatization process with 3 states.
% van der Pol ODEs for mu=1, with sample time 0.05s)
%
% Description:
% - Heptane to toluene aromatization process model described
%   by Watanabe and Himmelblau (1983).
%
% References:
% - Watanabe, K. and D. M. Himmelblau (1983). Fault diagnosis
%   in nonlinear chemical processes. AIChE J. 29, 243-261.
%     Part I. x(1)heory
%     Part II. Application to a Chemical Reactor
% - Robertson, D. G., and Lee, J. H. (1998). A Method for the
%   Estimation of Infrequent Abrupt Changes in Nonlinear
%   Systems.
%
% State variables
% xk(1) : x(1), reaction temperature, [K]
% xk(2) : Ch, outlet concentration of heptane, [gmol/m^3]
% xk(3) : Ct, outlet concentration of toluene, [gmol/m^3]
%
% Process inputs (disturbances)
% p(1) : k0, frequency factor, or pre-exponential factor [1/h]
% p(2) : U, overall heat transfer coefficient [J/(gmol.K)]
% 
% Output measurements
% y(1) : xk(1)
% y(2) : xk(3)
%
% Nominal values for the states and parameters
% x0 = [742.0;  % reaction temperature, [K]
%       463.0;  % outlet concentration of heptane, [gmol/m^3]
%       537.0];  % outlet concentration of toluene, [gmol/m^3]
% p0 = [5e8 / 1e8;  % (normalized)
%       6e5 / 1e5];  % (normalized)

    % Parameter values from normalized values in p(k)
    k0 = pk(1) .* 1e8;  % frequency or pre-exponential factor [h^-1]
    U = pk(2) .* 1e5;  % overall heat transfer coefficient [J/(gmol.K)]

    % Constant parameter values:
    Ea = params.Ea;  % activation energy (J/gmol)
    R = params.R;  % gas constant (J/(gmol.K))

    % Heat of reaction [J/gmol] as a function of temperature
    cD = params.cD;
    DeltaH = cD * [xk(1)^4; xk(1)^3; xk(1)^2; xk(1); 1];

    % Specific heat (molar)
    Cp = params.Cp;  % [J/(gmol.K)]

    % Other process parameters are assumed to be constant
    rho = params.rho;  % density, [gmol/m3]
    h = params.h;  % overall heat transfer coefficient, [J/(gmol.K)]
    A = params.A;  % area of heat exchange, [m^2]
    V = params.V;  % effective reactor volume, [m^3]
    q = params.q;  % inlet and outlet volumetric flow rate, [m^3/h]

    % Reaction rate constant
    %rate_const = k0 * exp(-Ea / (R * xk(1)));

    % Jacobian matrix
    dfdx = zeros(3, 3);

    dfdx(1, 1) = 1 - dt*((k0*xk(2)*exp(-Ea/(R*xk(1))) * ...
        (2*cD(3) + 6*cD(2)*xk(1) + 12*cD(1)*xk(1)^2))/(Cp*rho) ...
        + (Ea^2*k0*xk(2)*exp(-Ea/(R*xk(1))) * DeltaH) / (Cp*R^2*rho*xk(1)^4) ...
        + (2*Ea*k0*xk(2)*exp(-Ea/(R*xk(1))) * ...
            (cD(4) + 2*cD(3)*xk(1) + 4*cD(1)*xk(1)^3 + 3*cD(2)*xk(1)^2)) / ...
            (Cp*R*rho*xk(1)^2) ...
        - (2*Ea*k0*xk(2)*exp(-Ea/(R*xk(1))) * DeltaH) / (Cp*R*rho*xk(1)^3));

    dfdx(1, 2) = -dt*((k0*exp(-Ea/(R*xk(1))) * ...
            (cD(4) + 2*cD(3)*xk(1) + 4*cD(1)*xk(1)^3 + 3*cD(2)*xk(1)^2)) / (Cp*rho) ...
        + (Ea*k0*exp(-Ea/(R*xk(1)))*DeltaH) / (Cp*R*rho*xk(1)^2));

    dfdx(1, 3) = 0;

    dfdx(2, 1) = -(Ea*dt*k0*xk(2)*exp(-Ea/(R*xk(1))) * ...
            (Ea - 2*R*xk(1))) / (R^2*xk(1)^4);

    dfdx(2, 2) = 1 - (Ea*dt*k0*exp(-Ea/(R*xk(1))))/(R*xk(1)^2);

    dfdx(2, 3) = 0;

    dfdx(3, 1) = (Ea*dt*k0*xk(2)*exp(-Ea/(R*xk(1)))*(Ea - 2*R*xk(1))) / (R^2*xk(1)^4);

    dfdx(3, 2) = (Ea*dt*k0*exp(-Ea/(R*xk(1))))/(R*xk(1)^2);

    dfdx(3, 3) = 1;

end
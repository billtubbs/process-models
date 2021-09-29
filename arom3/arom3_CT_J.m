function J = arom3_CT_J(x, p, params)
% J = arom3_CT_J(x, p, params)
% Jacobian matrix of continous-time model of aromatization
% process with 3 states
%
% State variables
% x(1) : T, reaction temperature, [K]
% x(2) : Ch, outlet concentration of heptane, [gmol/m^3]
% x(3) : Ct, outlet concentration of toluene, [gmol/m^3]
%
% Process disturbances
% p(1) : k0, frequency factor, or pre-exponential factor [1/h]
% p(2) : U, overall heat transfer coefficient [J/(gmol.K)]

    % Parameter values from normalized values in p(k)
    k0 = p(1) .* 1e8;  % frequency or pre-exponential factor [h^-1]
    U = p(2) .* 1e5;  % overall heat transfer coefficient [J/(gmol.K)]

    % Constant parameter values:
    Ea = params.Ea;  % activation energy (J/gmol)
    R = params.R;  % gas constant (J/(gmol.K))

    % Heat of reaction [J/gmol] as a function of temperature
    %DeltaH = params.DeltaH(x(1));
    cD = params.cD;

    % Specific heat (molar)
    Cp = params.Cp;  % [J/(gmol.K)]

    % Other process parameters are assumed to be constant
    rho = params.rho;  % density, [gmol/m3]
    h = params.h;  % overall heat transfer coefficient, [J/(gmol.K)]
    A = params.A;  % area of heat exchange, [m^2]
    V = params.V;  % effective reactor volume, [m^3]
    q = params.q;  % inlet and outlet volumetric flow rate, [m^3/h]

    % Reaction rate constant
    rate_const = k0 * exp(-Ea / (R * x(1)));

%     % dT/dt, rate of change of reaction temperature
%     % (also x(1) in Watanbe and H.)
%     dx(1) = q / V * (u(1) - x(1)) - DeltaH / (rho * Cp) * rate_const * x(2) ...
%         + U * A / (rho * Cp * V) * (u(3) - x(1));
% 
%     % dCh/dt, rate of change of heptane concentration
%     % (x4 in Watanbe and H.)
%     dx(2) = q / V * (u(2) - x(2)) - rate_const * x(2);
% 
%     % dCt/dt, rate of change of toluene concentration
%     % (x(2) in Watanbe and H.)
%     dx(3) = -q / V * x(3) + rate_const * x(2);
% 

    % Jacobian matrix
    J = zeros(3, 3);

%     J(1, 1) = - q / V ...
%         - DeltaH * Ea * k0 * exp(-Ea/(R * x(1))) * x(2) / (Cp * rho * R * x(1)^2) ...
%         - U * A / (rho * Cp * V);
%     
%     J(1, 2) = - DeltaH / (rho * Cp) * rate_const;
%     
%     J(2, 1) = -Ea * k0 * exp(-Ea / (R * x(1))) / (R * x(1)^2) * x(2);
% 
%     J(2, 2) = -q / V - rate_const;
%     
%     J(3, 1) = Ea * k0 * exp(-Ea / (R * x(1))) / (R * x(1)^2) * x(2);
%     
%     J(3, 2) = rate_const;
% 
%     J(3, 3) = -q / V;

    % Recalculated with arom3_simplify.m:

    J(1, 1) = - q/V - (A*U)/(Cp*V*rho) ...
        - (k0*x(2)*exp(-Ea/(R*x(1)))*(cD(4) + 2*cD(3)*x(1) ...
        + 4*cD(1)*x(1)^3 + 3*cD(2)*x(1)^2))/(Cp*rho) ...
        - (Ea*k0*x(2)*exp(-Ea/(R*x(1))) * ...
            (cD(5) + cD(4)*x(1) + cD(1)*x(1)^4 + cD(2)*x(1)^3 + cD(3)*x(1)^2)) / ...
                (Cp*R*rho*x(1)^2);

    J(1, 2) = -(k0*exp(-Ea/(R*x(1))) * ...
        (cD(5) + cD(4)*x(1) + cD(1)*x(1)^4 + cD(2)*x(1)^3 + cD(3)*x(1)^2)) / (Cp*rho);

    J(1, 3) = 0;

    J(2, 1) = -(Ea*k0*x(2)*exp(-Ea/(R*x(1))))/(R*x(1)^2);

    J(2, 2) = - k0*exp(-Ea/(R*x(1))) - q/V;

    J(2, 3) = 0;

    J(3, 1) = (Ea*k0*x(2)*exp(-Ea/(R*x(1))))/(R*x(1)^2);

    J(3, 2) = k0*exp(-Ea/(R*x(1)));

    J(3, 3) = -q/V;

end
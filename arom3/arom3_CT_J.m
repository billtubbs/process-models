function J = arom3_CT_J(x, p, w, params)
% J = arom3_CT_J(x, p, w, params)
% Jacobian matrix of continous-time model of aromatization
% process with 3 states

    % State variables
    T = x(1);
    Ch = x(2);
    Ct = x(3);

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

    % Parameter values from normalized values in p(k)
    k0 = p(1) .* 1e8;  % frequency or pre-exponential factor [h^-1]
    U = p(2) .* 1e5;  % overall heat transfer coefficient [J/(gmol.K)]

    % Reaction rate constant
    rate_const = k0 * exp(-Ea / (R * x(1)));

%     % dT/dt, rate of change of reaction temperature
%     % (also x1 in Watanbe and H.)
%     dx(1) = q / V * (u(1) - x(1)) - DeltaH / (rho * Cp) * rate_const * x(2) ...
%         + U * A / (rho * Cp * V) * (u(3) - x(1));
% 
%     % dCh/dt, rate of change of heptane concentration
%     % (x4 in Watanbe and H.)
%     dx(2) = q / V * (u(2) - x(2)) - rate_const * x(2);
% 
%     % dCt/dt, rate of change of toluene concentration
%     % (x2 in Watanbe and H.)
%     dx(3) = -q / V * x(3) + rate_const * x(2);
% 

    % Jacobian matrix
    J = zeros(3, 3);

    J(1, 1) = - q / V ...
        - DeltaH * Ea * k0 * exp(-Ea/(R * x(1))) * x(2) / (Cp * rho * R * x(1)^2) ...
        - U * A / (rho * Cp * V);
    
    J(1, 2) = - DeltaH / (rho * Cp) * rate_const;
    
    J(2, 1) = -Ea * k0 * exp(-Ea / (R * x(1))) / (R * x(1)^2) * x(2);

    J(2, 2) = -q / V - rate_const;
    
    J(3, 1) = Ea * k0 * exp(-Ea / (R * x(1))) / (R * x(1)^2) * x(2);
    
    J(3, 2) = rate_const;

    J(3, 3) = -q / V;

end
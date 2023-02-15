function dydt = crane_ODEs(t,Y,params)
% dydt = crane_ODEs(t,Y,params)
% Right-hand-side of the ordinary differential equations (ODEs)
% for the dynamic model of a gantry crane system.
%
% This function can be used with MATLAB integration functions
% such as ode45 to simulate the system.
%
% Arguments:
%   t : double
%       Time.
%   Y : (6, 1) double vector
%       System outpus, y(t), defined:
%           y(1) : horizontal position of cart, x(t), left to right
%           y(2) : angle of pole, theta(t), clockwise from vertical
%           y(3) : downward force of cart on track, N_c(t)
%           y(4) : tension/compression in cable/arm, N(t)
%           y(5) : horizontal velocity of cart, Dxt(t), left to right
%           y(6) : angular velocity of pole, Dthetat(t), clockwise.
%   params : struct
%       Parameter values defined:
%           params.F : Force acting on cart (  +  ve to right)
%           params.L : Pendulum length
%           params.c_d : Drag coefficient of load
%           params.g : Acceleration due to gravity (  +  ve down)
%           params.m_c : Cart mass
%           params.m_p : Pole mass
%           params.muc : Coefficient of friction for cart and track 
%           params.r : Dadius of load
%           params.rho : Density of air.
% 
% The equations below were generated by the MATLAB
% Symbolic Math Toolbox version 9.0. 12 - Feb - 2023. See
% livescript 'solve_DAE_crane.mlx'.
% 

    % Output variables
    %x = Y(1,:);  % x(t) has no effect
    theta = Y(2,:);
    N_c = Y(3,:);
    %N = Y(4,:);  % this is recalculated
    Dxt = Y(5,:);
    Dthetat = Y(6,:);

    % Parameters
    F = params.F;
    L = params.L;
    c_d = params.c_d;
    g = params.g;
    m_c = params.m_c;
    m_p = params.m_p;
    muc = params.muc;
    r = params.r;
    rho = params.rho;

    % Temporary variables
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    mcmp = m_c + m_p;
    Dthetat_sq = Dthetat.^2;
    r_sq = r.^2;
    inv_mcmp = 1.0 ./ mcmp;

    % Components of velocity of load
    v_x = Dxt + Dthetat .* L .* cos_theta;
    v_y = Dthetat .* L .* sin_theta;
    v_x_sq = v_x.^2;
    v_y_sq = v_y.^2;

    % Components of aerodynamic drag force acting on load
    F_dx = (c_d .* rho .* r_sq .* v_x_sq .* pi .* sign(v_x)) ./ 2.0;
    F_dy = (c_d .* rho .* r_sq .* v_y_sq .* pi .* sign(v_y)) ./ 2.0;

    % Horizontal component of force on cable by cart
    %N_x = m_p .* ( ...
    %    D2xt2 + D2thetat2 .* L .* cos_theta - L .* sin_theta .* Dthetat_sq ...
    %);

    % Vertical component of force on cable by cart
    N_y = m_p .* ( ...
        g - L .* (D2thetat2 .* sin_theta + cos_theta .* Dthetat_sq) ...
    );

    % Downwards vertical force on track by cart
    N_c = N_y + g .* m_c;

    % Magnitude of tension (compression) in cable (arm)
    %N = sqrt(N_x .^ 2 + N_y .^ 2);

    % Friction force on cart by track
    F_f = N_c .* muc .* sign(Dxt .* N_c);

    % Clockwise angular acceleration of cable
    D2thetat2 = -(0.7500 * (sin_theta * (g * m_p - F_dy) - cos_theta^2 * ((m_p * (F - muc * sign(N_c(t) * Dxt) * N_c(t) + L * m_p * sin_theta * Dthetat_sq)) * inv_mcmp + F_dx))) / (L * m_p * ((0.7500 * m_p * cos_theta^2) * inv_mcmp - 1));

%     D2thetat2 = ( ...
%             ( ...
%                 cos_theta .* ( ...
%                     m_p .* inv_mcmp .* ( ...
%                         - F_f + F + L .* m_p .* sin_theta .* Dthetat_sq ...
%                     ) ...
%                     + F_dx ...
%                 ) ...
%                 - sin_theta .* (g .* m_p - F_dy) ...
%             ) .* (3.0 ./ 4.0) ...
%         ) ./ (L .* m_p .* (m_p .* cos_theta.^2 .* inv_mcmp .* (3.0 ./ 4.0) - 1.0));

    % Horizontal acceleration of cart
    D2xt2 = - inv_mcmp .* ( ...
        F ...
        - F_f ...
        - L .* m_p .* (D2thetat2 .* cos_theta - sin_theta .* Dthetat_sq) ...
    );

    % Return vector of time derivatives of y(t)
    dydt = [
        Dxt;
        Dthetat;
        N_c;
        D2xt2;
        D2thetat2;
    ];

end
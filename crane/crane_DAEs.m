function dy = crane_DAEs(t,Y,YP,params)
% dy = crane_DAEs(t,Y,YP,params)
% Differential algebraic equations (DAEs) for the 
% dynamic model of a gantry crane system.
%
% This function can be used with the MATLAB ode15i
% solver to simulate the system.
%
% Arguments:
%   t : double
%       Time.
%   Y : (5, 1) double vector
%       System outpus, y(t), defined:
%           y(1) : horizontal position of cart, x(t), from left to right
%           y(2) : angle of pole, theta(t), clockwise from vertical
%           y(3) : downward force of cart on track, Nc(t)
%           y(4) : velocity of cart, Dxt(t), from left to right
%           y(5) : angular velocity of pole, Dthetat(t), clockwise.
%   YP : (5, 1) double vector
%       Partial derivatives of y(t).
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
    %x = Y(1,:);  % x(t) is not used
    %theta = Y(2,:);
    %Nc = Y(3,:);
    %Dxt = Y(4,:);
    %Dthetat = Y(5,:);

    % Partial derivatives
    %YP1 = YP(1,:);
    %YP2 = YP(2,:);
    % YP3 = YP(3,:)  % not used
    %YP4 = YP(4,:);
    %YP5 = YP(5,:);

    F_d = Y(5,:);
    F_f = Y(4,:);
    N_c = Y(3,:);
    N_x = Y(6,:);
    N_y = Y(7,:);
    v_x = Y(8,:);
    v_y = Y(9,:);
    Dxt = Y(10,:);
    Dthetat = Y(11,:);
    
    YP1 = YP(1,:);
    YP2 = YP(2,:);
    YP10 = YP(10,:);
    YP11 = YP(11,:);
    theta = Y(2,:);

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

    cos_theta = cos(theta);
    sin_theta = sin(theta);
    Dthetat_sq = Dthetat.^2;
    r_sq = r.^2;
    v_x_sq = v_x.^2;
    v_y_sq = v_y.^2;

    % Compute DAEs
    dy = nan(11, 1);

    dy(1) = - Dxt + v_x - Dthetat .* L .* cos_theta;

    dy(2) = v_y - Dthetat .* L .* sin_theta;

    dy(3) = F_d - (c_d .* rho .* r_sq .* pi .* (v_x_sq + v_y_sq)) ./ 2.0;

    dy(4) = N_x - m_p .* ( ...
        YP10 + YP11 .* L .* cos_theta - L .* sin_theta .* Dthetat_sq ...
    );

    dy(5) = N_y - m_p .* ( ...
        g - L .* (YP11 .* sin_theta + cos_theta .* Dthetat_sq) ...
    );

    dy(6) = N_c - N_y - g .* m_c;

    dy(7) = F_f - N_c .* muc .* sign(Dxt .* N_c);

    dy(8) = YP10 + ( ...
        F_f ...
        - F ...
        + L .* m_p .* (YP11 .* cos_theta - sin_theta .* Dthetat_sq) ...
    ) ./ (m_c + m_p);

    dy(9) = YP11 + ( ...
            ( ...
            cos_theta .* ( ...
                YP10 .* m_p ...
                + (c_d .* rho .* r_sq .* v_x_sq .* pi .* sign(v_x)) ./ 2.0 ...
            ) ...
            - sin_theta .* ( ...
                g .* m_p ...
                - (c_d .* rho .* r_sq .* v_y_sq .* pi .* sign(v_y)) ./ 2.0 ...
            ) ...
        ) .* (3.0 ./ 4.0) ...
    ) ./ (L .* m_p);

    dy(10) = Dxt - YP1;

    dy(11) = Dthetat - YP2;

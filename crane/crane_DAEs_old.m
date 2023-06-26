function dy = crane_DAEs_old(t,Y,YP,params)
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
%           params.F : Force acting on cart ( + ve to right)
%           params.L : Pendulum length
%           params.c_d : Drag coefficient of load
%           params.g : Acceleration due to gravity ( + ve down)
%           params.m_c : Cart mass
%           params.m_p : Pole mass
%           params.muc : Coefficient of friction for cart and track 
%           params.r : Dadius of load
%           params.rho : Density of air.
% 
% The equations below were generated by the MATLAB
% Symbolic Math Toolbox version 9.0. 25 - Jan - 2023. See
% livescript 'solve_DAE_crane.mlx'.
% 

    % Output variables
    %x = Y(1,:);  % x(t) is not used
    theta = Y(2,:);
    Nc = Y(3,:);
    Dxt = Y(4,:);
    Dthetat = Y(5,:);

    % Partial derivatives
    YP1 = YP(1,:);
    YP2 = YP(2,:);
    % YP3 = YP(3,:)  % not used
    YP4 = YP(4,:);
    YP5 = YP(5,:);

    % Parameters
    F = params.F;
    L = params.L;
    cd = params.c_d;
    g = params.g;
    mc = params.m_c;
    mp = params.m_p;
    muc = params.muc;
    r = params.r;
    rho = params.rho;
    
    % Temporary variables
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    mcmp = mc + mp;
    Dthetat_sq = Dthetat.^2;
    DxtNc = Dxt .* Nc;
    sgn_DxtNc = sign(DxtNc);
    inv_mcmp = 1.0 ./ mcmp;
    pi3div2 = 1.550313834014991e+1;

    % Compute DAEs
    dy = nan(5, 1);

    dy(1) = Nc .* 1.0e+3 - g .* mcmp ...
        + L .* mp .* (YP5 .* sin_theta + cos_theta .* Dthetat_sq);

    dy(2) = YP5 - ( ...
            -g .* sin_theta ...
            + cos_theta .* ( ...
                inv_mcmp .* ( ...
                    F + ...
                    L .* mp .* Dthetat_sq .* ( ...
                        sin_theta + muc .* sgn_DxtNc .* cos_theta ...
                    ) ...
                ) ...
                - g .* muc .* sgn_DxtNc ...
            ) ...
            + ( ...
                L.^2 .* cd .* r.^2 .* rho .* Dthetat_sq .* pi3div2 ...
            ) ./ mp ...
        ) ./ ( ...
            L .* ( ...
                mp .* cos_theta .* inv_mcmp .* ( ...
                    cos_theta - muc .* sgn_DxtNc ...
                ) ...
                - 4.0 ./ 3.0 ...
            ) ...
        );

    dy(3) = YP4 + inv_mcmp .* ( ...
            -F ...
            + L .* mp .* (YP5 .* cos_theta - sin_theta .* Dthetat_sq) ...
            + Nc .* muc .* sgn_DxtNc .* 1.0e+3 ...
        );

    dy(4) = Dxt - YP1;

    dy(5) = Dthetat - YP2;

end
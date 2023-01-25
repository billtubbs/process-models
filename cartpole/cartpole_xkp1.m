function xkp1 = cartpole_xkp1(xk,uk,params)
% xkp1 = cartpole_xkp1(xk,uk,params)
% State transition function for a discrete-time 
% model of the non-linear cart-pole system from
% the DAEs defined by Florian (2007).
%
% See script 'cartpole_DAEs.m' for source of eqn's below.
%
% Arguments:
%   xk : (4, 1) double vector
%       System states at time k, defined:
%           xk(1) : position of cart, x(t), from left to right
%           xk(2) : velocity of cart, Dx(t), from left to right
%           xk(3) : angle of pole, theta(t), clockwise from vertical
%           xk(4) : angular velocity of pole, Dtheta(t), clockwise.
%   uk : (1, 1) double
%       Input variable : force acting on cart (+ve to right).
%   params : struct
%       Parameter values defined:
%           params.g : acceleration due to gravity (+ve down)
%           params.l : pendulum length
%           params.mc : cart mass
%           params.mp : pole mass
%           params.muc : coefficient of friction for cart and track 
%           params.mup : coefficient of friction for pole joint
%           params.dt : sampling interval.
%

    assert(isequal(size(xk),[4 1]))
    assert(isscalar(uk))

    % Sampling interval
    dt = params.dt;

    % System variables at time k
    %xt = xk(1);
    Dxt = xk(2);
    thetat = xk(3);
    Dthetat = xk(4);

    % Set force on cart, F, to value of input variable
    F = uk;

    % Get parameter values
    g = params.g;
    l = params.l;
    mc = params.mc;
    mp = params.mp;
    muc = params.muc;
    mup = params.mup;

    % Temporary variables
    cos_theta = cos(thetat);
    sin_theta = sin(thetat);
    mcmp = mc + mp;
    Dthetat_sq = Dthetat.^2;
    invl = 1.0 ./ l;
    %sgn_DxtNc = sign(Dxt .* Nc);  % see below
    inv_mcmp = 1.0 ./ mcmp;

    % Compute Nc and derivatives of system variables
    % Nc is defined as the downward force on the track
    % by the cart and is used to compute the friction.

    % First, compute second derivative of theta (i.e. 
    % angular acceleration) assuming Nc is positive, where
    % Nc is the downward force of cart on track.
    % Usually, for common choices of the parameters, Nc 
    % will be always positive (Florian, 2007).
    sgn_DxtNc = sign(Dxt);
    D2thetat2 = calculate_D2thetat2(Dthetat, Dthetat_sq, F, g, l, mp, ...
        muc, mup, sin_theta, cos_theta, invl, inv_mcmp, sgn_DxtNc);

    % Now compute Nc using this D2thetat value
    Nc = (g .* mcmp - l .* mp ...
        .* (D2thetat2 .* sin_theta + Dthetat_sq .* cos_theta) ...
    );

    % If necessary, re-compute D2thetat2
    %assert(sign(Nc) == 1)
    if sign(Nc) == -1
        % TODO: The only way to avoid this would be to 
        % save Nc as an additional state.
        sgn_DxtNc = -sgn_DxtNc;
        D2thetat2 = calculate_D2thetat2(Dthetat, Dthetat_sq, F, g, l, ...
            mp, muc, mup, sin_theta, cos_theta, invl, inv_mcmp, sgn_DxtNc);
    end

    D2xt2 = inv_mcmp .* ( ...
            F + mp .* l .* (Dthetat_sq .* sin_theta - D2thetat2 .* cos_theta) ...
            - muc .* Nc .* sgn_DxtNc ...
        );

    % Partial derivatives at time t
    dy = [
        Dxt;
        D2xt2;
        Dthetat;
        D2thetat2
    ];

    % x(k+1) by Euler method
    xkp1 = xk + dt * dy;

    % Make sure theta stays within -pi and pi
    xkp1(3) = wrapToPi(xkp1(3));

end


function D2thetat = calculate_D2thetat2(Dthetat, Dthetat_sq, F, g, l, mp, ...
    muc, mup, sin_theta, cos_theta, invl, inv_mcmp, sgn_DxtNc)
% Compute second derivative of pole angle, D^2*theta(t).

    D2thetat = ( ...
            invl .* ( ...
                g .* sin_theta ...
                + cos_theta .* ( ...
                    inv_mcmp .* ( ...
                        -F ...
                        - mp .* l .* Dthetat_sq .* ( ...
                            sin_theta ...
                            + muc .* cos_theta .* sgn_DxtNc ...
                        ) ...
                    ) ...
                    + muc .* g .* sgn_DxtNc ...
                ) ...
            - (mup .* Dthetat .* invl) ./ mp) ...
        ) ./ ( ...
            4.0 ./ 3.0 ...
            - mp .* cos_theta .* inv_mcmp .* (cos_theta - muc .* sgn_DxtNc) ...
        );

end
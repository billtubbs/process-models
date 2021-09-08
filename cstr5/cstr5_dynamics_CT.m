%% Continuous-time nonlinear model of CSTR reactor

function dx = cstr5_dynamics_CT(x, p, w, params)
    % CSTR dynamic model
    %
    % Description:
    % - Constinuously-stirred reactor (CSTR) system for the free-radical
    %   polymerization of methyl methacrylate (MMA) with azo-bis-
    %   isobutyronitrile (AIBN) as initiator and toluene as solvent.
    %
    % See :
    % - Daoutidis, Soroush, and Kravaris, 1990, and
    % - Robertson, Kesavan and Lee, 1995
    %
    % State variables
    % x(1) : Cm
    % x(2) : CI
    % x(3) : T
    % x(4) : D0
    % x(5) : DI
    %
    % Process disturbances
    % p(1) : CI_in
    % p(2) : CM_in
    % 
    % Nominal values for the states and parameters
    % (i) as published in the paper:
    % x0 = [5.53;  % kmol.m^-3
    %       0.684;  % kmol.m^-3
    %       331.8;  % K
    %       0.0019;  % kg.m^-3
    %       47.4  % kg.m^-3
    %       ];
    % (ii) based on this implementation:
    % x0 = [5.630;  % kmol.m^-3
    %       0.6388;  % kmol.m^-3
    %       331.3;  % K
    %       0.0017;  % kg.m^-3
    %       44.27  % kg.m^-3
    %       ];
    % p0 = [8.0;  % kmol.m^-3
    %       6.6;  % kmol.m^-3
    %       ];

    % Kinetic parameters
    Z.Tc = params.Z.Tc;
    Z.Td = params.Z.Td;  % kmol/m^3.h
    Z.I = params.Z.I;  % h^-1
    Z.P = params.Z.P;  % kmol/m^3.h
    Z.fm = params.Z.fm;  % kmol/m^3.h
    E.Tc = params.E.Tc;  % kJ/kmol
    E.Td = params.E.Td;  % kJ/kmol
    E.I = params.E.I;  % kJ/kmol
    E.P = params.E.P;  % kJ/kmol
    E.fm = params.E.fm;  % kJ/kmol
    
    % Physical parameters
    F = params.F;  % m^3.h^-1
    V = params.V;  % m^3
    rho = params.rho;  % kg.m^-3
    Tj = params.Tj;  % K
    Tin = params.Tin;  % K
    A = params.A;  % m^2
    f_star = params.f_star;
    FI = params.FI;  % m^3.h^-1
    Fm = params.Fm;  % m^3.h^-1
    Mm = params.Mm;  % kg.kmol^-1
    mDeltaHp = params.mDeltaHp;  % kJ.kmol^-1
    cp =  params.cp;  % kJ.kg^-1
    U = params.U;  % kJ.h^-1.K-1.m^-2
    R = params.R;  % kJ.mol^-1.K^-1
    
    % Disturbances
    Ci_in = p(1);  % inlet initiator concentration
    Cm_in = p(2);  % inlet monomer concentration
    
    % Intermediate calculations
    Zp_exp_Ep = Z.P * exp( -E.P / (R*x(3)) );
    Zfm_exp_Efm = Z.fm * exp( -E.fm / (R*x(3)) );
    ZTd_exp_ETd = Z.Td * exp( -E.Td / (R*x(3)) );
    ZTc_exp_ETc = Z.Tc * exp( -E.Tc / (R*x(3)) );
    P0_CI = sqrt( ...
        2 * f_star * x(2) * Z.I * exp( -E.I / (R*x(3)) ) ...
        / ( ZTd_exp_ETd + ZTc_exp_ETc ) ...
    );

    % Add process disturbances to each state
    x = x + w;

    % State equations
    dx = nan(5, 1);

    % dCm/dt
    dx(1) = -( Zp_exp_Ep + Zfm_exp_Efm ) * x(1) * P0_CI ...
        + (Fm * Cm_in - F * x(1)) / V;

    % dCI/dt
    dx(2) = -Z.I * exp( -E.I / (R*x(3)) ) * x(2) ...
        + (FI * Ci_in - F * x(2)) / V;

    % dT/dt
    dx(3) = Zp_exp_Ep * x(1) * mDeltaHp / (rho * cp) * P0_CI ...
        - U * A / (rho * cp * V) * (x(3) - Tj) + F*(Tin - x(3)) / V;

    % dD0/dt
    dx(4) = ( ZTc_exp_ETc / 2 + ZTd_exp_ETd ) * P0_CI^2 ...
        + Zfm_exp_Efm * x(1) * P0_CI - F * x(4) / V;

    % dD1/dt
    dx(5) = Mm * ( Zp_exp_Ep + Zfm_exp_Efm ) * x(1) * P0_CI ...
        - F * x(5) / V;

end
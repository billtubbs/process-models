%% Discrete-time state equations for nonlinear model of CSTR reactor

function x = cstr5_dynamics(x, p, w, params)
    % CSTR model
    %
    % See :
    % - Daoutidis, Soroush, and Kravaris, 1990, and
    % - Robertson, Kesavan and Lee, 1995
    % 
    % Description:
    % - Constinuously-stirred reactor (CSTR) system for the free-radical
    %   polymerization of methyl methacrylate (MMA) with azo-bis-
    %   isobutyronitrile (AIBN) as initiator and toluene as solvent.
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

    % Continuous-time differential equations
    dx = cstr5_dynamics_CT(x, p, w, params);

    % Sample time (hours)
    dt = 0.016667; 
    
    % Euler integration of continuous-time dynamics x' = f(x)
    x = x + dx * dt;

end
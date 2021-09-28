function [dx, y] = cstr5(t, x, p, w, params, varargin)
% [dx, y] = cstr5(t, x, p, w, params, varargin)
% Continous-time model of CSTR process with 5 states
%
% Description:
% - Constinuously-stirred reactor (CSTR) system for the free-radical
%   polymerization of methyl methacrylate (MMA) with azo-bis-
%   isobutyronitrile (AIBN) as initiator and toluene as solvent.
%
% References:
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
%
% This function follows the conventions of grey-box model
% files in MATLAB except that the parameters have been grouped
% into a struct (therefore can't be estimated).
%
% See documentation on defining grey-box model files:
% https://www.mathworks.com/help/ident/ug/creating-idnlgrey-model-files.html

    % RHS of differential equation
    dx = cstr5_dynamics_CT(x, p, w, params);
    
    % Output equations
    y = cstr5_outputs(x, p);

end
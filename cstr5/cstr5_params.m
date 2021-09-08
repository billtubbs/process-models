% Constants and parameters for CSTR reactor simulation model
%
% See :
% - Daoutidis, Soroush, and Kravaris, 1990, and
% - Robertson, Kesavan and Lee, 1995

% Kinetic parameters
params.Z.Tc = 3.8223e10;  % kmol/m^3.h
params.Z.Td = 3.1457e11;  % kmol/m^3.h
params.Z.I = 3.7920e18;  % h^-1
params.Z.P = 1.7700e9;  % kmol/m^3.h
params.Z.fm = 1.0067e15;  % kmol/m^3.h
params.E.Tc = 2.9442e3;  % kJ/kmol
params.E.Td = 2.9442e3;  % kJ/kmol
params.E.I = 1.2877e5;  % kJ/kmol
params.E.P = 1.8283e4;  % kJ/kmol
params.E.fm = 7.4478e4;  % kJ/kmol

% Physical parameters
params.F = 1.0;  % m^3.h^-1
params.V = 0.1;  % m^3
params.rho = 866;  % kg.m^-3
params.Tj = 295;  % K
params.Tin = 350;  % K
params.A = 2.0;  % m^2
params.f_star = 0.58;
params.FI = 0.08;  % m^3.h^-1
params.Fm = 0.92;  % m^3.h^-1
params.Mm = 100.12;  % kg.kmol^-1
params.mDeltaHp = 57800;  % kJ.kmol^-1
params.cp =  2.0;  % kJ.kg^-1
params.U = 800;  % kJ.h^-1.K-1.m^-2
params.R = 8.314;  % kJ.mol^-1.K^-1

% Nominal values for the states and parameters
% (i) as published in the paper:
% x0 = [5.53;   % molar concentration of monomer (kmol.m^-3)
%       0.684;  % molar concentration of initiator (kmol.m^-3)
%       331.8;  % reactor temperature (K)
%       0.0019;  % molar concentration of dead polymer chains (kmol.m^-3)
%       47.4   % mass concentration of dead polymer chains (kg.m^-3)
%       ];
% (ii) Steady-state based on numerically solving eqns:
x0 = [5.63047;   % molar concentration of monomer (kmol.m^-3)
      0.63880;   % molar concentration of initiator (kmol.m^-3)
      331.2426;  % reactor temperature (K)
      0.0016561;  % molar concentration of dead polymer chains (kmol.m^-3)
      44.2056;  % mass concentration of dead polymer chains (kg.m^-3)
      ];
p0 = [8.0;  % molar concentration of monomer in the monomer inlet stream (kmol.m^-3)
      6.6;  % molar concentration of initiator in the initiator inlet stream (kmol.m^-3)
      ];

% Deriving linear approximations of CSTR model symbolically
clear

% Kinetic parameters
syms Z_Tc Z_Td Z_I Z_P Z_fm E_Tc E_Td E_I E_P E_fm

% Physical parameters
syms F V rho Tj Tin A f_star FI Fm Mm mDeltaHp cp U R

% State variables including disturbance parameters
syms Cm CI T D0 D1 Ci_in Cm_in
x = [Cm;
     CI;
     T;
     D0;
     D1;
     Ci_in;
     Cm_in];

% Intermediate calculations
Zp_exp_Ep = Z_P * exp( -E_P / (R*x(3)) );
Zfm_exp_Efm = Z_fm * exp( -E_fm / (R*x(3)) );
ZTd_exp_ETd = Z_Td * exp( -E_Td / (R*x(3)) );
ZTc_exp_ETc = Z_Tc * exp( -E_Tc / (R*x(3)) );
ZI_expEI = Z_I * exp( -E_I / (R*x(3)) );
P0_CI = sqrt( ...
    2 * f_star * x(2) * Z_I * exp( -E_I / (R*x(3)) ) ...
    / ( ZTd_exp_ETd + ZTc_exp_ETc ) ...
);

% State equations
% dx(1)
dCm = -( Zp_exp_Ep + Zfm_exp_Efm ) * x(1) * P0_CI ...
    + (Fm * Cm_in - F * x(1)) / V;

% dx(2)
dCI = -ZI_expEI * x(2) ...
    + (FI * Ci_in - F * x(2)) / V;

% dx(3)
dT = Zp_exp_Ep * x(1) * mDeltaHp / (rho * cp) * P0_CI ...
    - U * A / (rho * cp * V) * (x(3) - Tj) + F*(Tin - x(3)) / V;

% dx(4)
dD0 = ( ZTc_exp_ETc / 2 + ZTd_exp_ETd ) * P0_CI^2 ...
    + Zfm_exp_Efm * x(1) * P0_CI - F * x(4) / V;

% dx(5)
dD1 = Mm * ( Zp_exp_Ep + Zfm_exp_Efm ) * x(1) * P0_CI ...
    - F * x(5) / V;

% dx(6)
dCi_in = 0;

% dx(7)
dCm_in = 0;

dx = [dCm;
      dCI;
      dT;
      dD0;
      dD1;
      dCi_in;
      dCm_in];
n = size(x,1);
assert(size(dx,1) == n)

% Compute Jacobian
dfdx = sym('dfdx', [n n]);
for i=1:n
    for j=1:n
        dfdx(i,j) = simplify(diff(dx(i), x(j)));
        fprintf("\ndfdx(%d, %d) = ", i, j)
        disp(dfdx(i,j))
    end
end
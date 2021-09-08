function J = cstr7_CT_J(x, w, params)
    % Jacobian matrix of continous-time CSTR model based on :
    % - Daoutidis, Soroush, and Kravaris, 1990 and
    %   presented by Robertson, Kesavan and Lee, 1995

    % State variables
    Cm = x(1);
    CI = x(2);
    T = x(3);
    %D0 = x(4);  % not needed
    %D1 = x(5);  % not needed
    %Ci_in = x(6);  % not needed
    %Cm_in = x(7);  % not needed

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
    %Tj = params.Tj;  % K
    %Tin = params.Tin;  % K
    A = params.A;  % m^2
    f_star = params.f_star;
    FI = params.FI;  % m^3.h^-1
    Fm = params.Fm;  % m^3.h^-1
    Mm = params.Mm;  % kg.kmol^-1
    mDeltaHp = params.mDeltaHp;  % kJ.kmol^-1
    cp =  params.cp;  % kJ.kg^-1
    U = params.U;  % kJ.h^-1.K-1.m^-2
    R = params.R;  % kJ.mol^-1.K^-1

    % Intermediate calculations
    Zp_exp_Ep = Z.P * exp(-E.P / (R*T));
    Zfm_exp_Efm = Z.fm * exp(-E.fm / (R*T));
    ZTd_exp_ETd = Z.Td * exp(-E.Td / (R*T));
    ZTc_exp_ETc = Z.Tc * exp(-E.Tc / (R*T));
    ZI_expEI = Z.I * exp(-E.I / (R*T));
    P0_CI = sqrt((2 * CI * f_star * ZI_expEI) / (ZTc_exp_ETc + ZTd_exp_ETd));

    % Jacobian matrix
    J = zeros(7, 7);

    J(1, 1) = -(Zp_exp_Ep + Zfm_exp_Efm) * P0_CI - F/V;

    J(1, 2) = -(Cm * f_star * ZI_expEI * (Zp_exp_Ep + Zfm_exp_Efm)) ...
        / ((ZTc_exp_ETc + ZTd_exp_ETd) * P0_CI);

    J(1, 3) = - (Cm * P0_CI * (E.fm * Zfm_exp_Efm + E.P * Zp_exp_Ep)) / (R*T^2) ... 
        - (CI * Cm * ZI_expEI * f_star * (Zfm_exp_Efm + Zp_exp_Ep) ...
           * (E.I * ZTc_exp_ETc + E.I * ZTd_exp_ETd - E.Tc * ZTc_exp_ETc - E.Td * ZTd_exp_ETd)) ... 
        / (P0_CI * R*T^2 * (ZTc_exp_ETc + ZTd_exp_ETd)^2);
    
    J(1, 7) = Fm/V;

    J(2, 2) = -ZI_expEI - F/V;

    J(2, 3) = -(CI * E.I * ZI_expEI) / (R*T^2);
    
    J(2, 6) = FI/V;

    J(3, 1) = (mDeltaHp * Zp_exp_Ep * P0_CI) / (cp * rho);

    J(3, 2) = (Cm * f_star * mDeltaHp * ZI_expEI * Zp_exp_Ep) ...
        / (cp * rho * (ZTc_exp_ETc + ZTd_exp_ETd) * P0_CI);

    J(3, 3) = (Cm * E.P * P0_CI * Zp_exp_Ep * mDeltaHp) ...
        / (R*T^2 * cp * rho) ...
        - (A * U) / (V * cp * rho) ...
        - F/V ...
        + (CI * Cm * ZI_expEI * Zp_exp_Ep * f_star * mDeltaHp ...
           * (E.I * ZTc_exp_ETc + E.I * ZTd_exp_ETd - E.Tc * ZTc_exp_ETc - E.Td * ZTd_exp_ETd)) ...
        / (P0_CI * R*T^2 * cp * rho * (ZTc_exp_ETc + ZTd_exp_ETd)^2);

    J(4, 1) = Zfm_exp_Efm * P0_CI;

    J(4, 2) = (ZI_expEI * f_star * (P0_CI * ZTc_exp_ETc + 2 * P0_CI*ZTd_exp_ETd + Cm * Zfm_exp_Efm)) ...
        / (P0_CI * (ZTc_exp_ETc + ZTd_exp_ETd));

    J(4, 3) = (Cm * E.fm * P0_CI * Zfm_exp_Efm) / (R*T^2) ...
        + (CI * E.I * ZI_expEI * f_star * (ZTc_exp_ETc + 2 * ZTd_exp_ETd)) ...
        / (R*T^2 * (ZTc_exp_ETc + ZTd_exp_ETd)) ...
        - (CI * ZI_expEI * f_star * (E.Tc * ZTc_exp_ETc + E.Td * ZTd_exp_ETd) ...
           * (ZTc_exp_ETc + 2 * ZTd_exp_ETd)) ...
        / (R*T^2 * (ZTc_exp_ETc + ZTd_exp_ETd)^2) ...
        + (CI * ZI_expEI * Z.Td * f_star * (E.Tc * ZTc_exp_ETc + 2 * E.Td * ZTd_exp_ETd)) ...
        / (R*T^2 * (ZTd_exp_ETd + ZTc_exp_ETc * Z.Td)) ...
        + (CI * Cm * ZI_expEI * Zfm_exp_Efm * f_star ...
           * (E.I * ZTc_exp_ETc + E.I * ZTd_exp_ETd - E.Tc * ZTc_exp_ETc - E.Td * ZTd_exp_ETd)) ...
        / (P0_CI * R*T^2 * (ZTc_exp_ETc + ZTd_exp_ETd)^2);

    J(4, 4) = -F/V;

    J(5, 1) = Mm * (Zp_exp_Ep + Zfm_exp_Efm) * P0_CI;

    J(5, 2) = (Cm * Mm * f_star * ZI_expEI * (Zp_exp_Ep + Zfm_exp_Efm)) ...
        / ((ZTc_exp_ETc + ZTd_exp_ETd) * P0_CI);

    J(5, 3) = (Cm * Mm * P0_CI * (E.fm * Zfm_exp_Efm + E.P * Zp_exp_Ep)) ...
        / (R*T^2) ...
        + (CI * Cm * Mm * ZI_expEI * f_star * (Zfm_exp_Efm + Zp_exp_Ep) ...
        * (E.I * ZTc_exp_ETc + E.I * ZTd_exp_ETd - E.Tc * ZTc_exp_ETc - E.Td * ZTd_exp_ETd)) ...
        / (P0_CI * R*T^2 * (ZTc_exp_ETc + ZTd_exp_ETd)^2);

    J(5, 5) = -F/V;

end
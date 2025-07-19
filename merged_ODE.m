%merged_ODE.m
function dydt = merged_ODE(t, y, t_f, app)
% tableA5_kinetics - Kinetic model for Chen's Table A-5 reactions (1 to 70)
%
% This function computes the time-derivatives for the 70 reactions in
% Chen’s Table A-5 (combined calibrated breakpoint, UV/HOCl, and UV/chloramine).

%
% Photolysis-rates are based on the formula:
%   r = (E0/L)*(1-10^(-a_total*L))*(a_species/a_total)*phi
% with :
% E0: irradiance (mol/m2.s)
% a: absorbances (M⁻¹cm⁻¹):
% L: optical pathlength (m)
%
% r = (E0/L)*(1-10^(-a_total*L))*(a_species/a_total)*phi
% with:

% All reaction rates and stoichiometries are taken directly from Table A-5.
%
% IMPORTANT MODIFICATION: The reaction “•OH + OH⁻ → O•⁻ + H2O” (Reaction 40)
% and the corresponding consumption of O•⁻ in Table C-1 (Reaction 16) are now
% implemented so that the untracked O•⁻ is now tracked as species 49 (O_rad_minus).
%
% The final derivative vector is multiplied by t_f.(app.t_unit).

%% Extract species from state vector (61×1) 
% all species should be in mol/L at this point
NH2Cl         = y(1);
NHCl2         = y(2);
HOCl          = y(3);
OCl_minus     = y(4);
NHCl_rad      = y(5);
NH2_rad       = y(6);
NCl2_rad      = y(7);
NCl3          = y(8);
NH3           = y(9);
NH4           = y(10);
H2O2          = y(11);
OH_rad        = y(12);
HO2_rad       = y(13);
HO2_minus     = y(14);
O3P           = y(15);
O1D           = y(16);
Cl_rad        = y(17);
Cl2_rad_minus = y(18);
Cl2           = y(19);
ClO_rad       = y(20);
ClO2_minus    = y(21);
ClO2_rad      = y(22);
ClO3_minus    = y(23);
Cl2O2         = y(24);
ClOH_rad_minus = y(25);
Cl3_minus     = y(26);
NOH           = y(27);
NO_rad        = y(28);
NO2_rad       = y(29);
NO2_minus     = y(30);
N2O3          = y(31);
N2O4          = y(32);
NO3_minus     = y(33);
dioxane       = y(34);
DOC           = y(35);
O2            = y(36);
H_plus        = y(37);
OH_minus      = y(38);
HCO3_minus    = y(39);
H2CO3         = y(40);
CO3_2         = y(41);
Cl_minus      = y(42);
HCl           = y(43);
NH2O2_rad     = y(44);
NHClO2_rad    = y(45);
N2O           = y(46);
CO3_rad       = y(47);
N2            = y(48);
O_rad_minus   = y(49);  % Newly tracked O•⁻
HPO4_2        = y(50);  % HPO4²⁻
H2PO4_minus   = y(51);  % H2PO4⁻
HPO4_rad      = y(52);  % phosphate radical HPO4•⁻
HClOH_rad     = y(53); 
O2_super      = y(54); 
O3_var =y(55);          %MOLECULAR ozone
O3_rad = y(56);          %O3 radical
H_atom = y(57); %H atom
ClO2 = y(58); %ClO2
O3_rad_minus = y(59); % O3•⁻
NOM = y(60); %NOM
caffeine = y(61) ; %caffeine
sucralose = y(62);


%% Initialize derivative variables (all set to zero)
dNH2Cl_dt         = 0;
dNHCl2_dt         = 0;
dHOCl_dt          = 0;
dOCl_minus_dt     = 0;
dNHCl_rad_dt      = 0;
dNH2_rad_dt       = 0;
dNCl2_rad_dt      = 0;
dNCl3_dt          = 0;
dNH3_dt           = 0;
dNH4_dt           = 0;
dH2O2_dt          = 0;
dOH_rad_dt        = 0;
dHO2_rad_dt       = 0;
dHO2_minus_dt     = 0;
dO3P_dt           = 0;
dO1D_dt           = 0;
dCl_rad_dt        = 0;
dCl2_rad_minus_dt = 0;
dCl2_dt           = 0;
dClO_rad_dt       = 0;
dClO2_minus_dt    = 0;
dClO2_rad_dt      = 0;
dClO3_minus_dt    = 0;
dCl2O2_dt         = 0;
dClOH_rad_minus_dt = 0;
dCl3_minus_dt     = 0;
dNOH_dt           = 0;
dNO_rad_dt        = 0;
dNO2_rad_dt       = 0;
dNO2_minus_dt     = 0;
dN2O3_dt          = 0;
dN2O4_dt          = 0;
dNO3_minus_dt     = 0;
ddioxane_dt       = 0;
dDOC_dt           = 0;
dO2_dt            = 0;
dH_plus_dt        = 0;
dOH_minus_dt      = 0;
dHCO3_minus_dt    = 0;
dH2CO3_dt         = 0;
dCO3_2_dt         = 0;
dCl_minus_dt      = 0;
dHCl_dt           = 0;
dNH2O2_rad_dt     = 0;
dNHClO2_rad_dt    = 0;
dN2O_dt           = 0;
dCO3_rad_dt       = 0;
dN2_dt            = 0;
dO_rad_minus_dt   = 0;  % For species 49
dHPO4_2_dt        = 0;  % For species 50
dH2PO4_minus_dt   = 0;  % For species 51
dHPO4_rad_dt      = 0;  % For species 52
dHClOH_rad_dt = 0;
dO2_super_dt     = 0;
dO3_var_dt =0;
dO3_rad_dt = 0;
dH_atom_dt = 0; %H atom
dClO2_dt = 0; %ClO2
dO3_rad_minus_dt = 0; %O3 radical
dNOM_dt = 0; %NOM
dcaffeine_dt = 0; %caffeine
dsucralose_dt = 0;
%% Global parameters and photolysis definitions
T = app.T;    % Temperature (K)
L = app.L;        % Fixed optical path length (cm)
V = app.V;   % Illuminated volume (L)
A = app.A;   % Illuminated Area (m^2)
E0 = app.E0;  % incident light intensity in mol photons/cm2/s

% Absorbance coefficients (M⁻¹cm⁻¹):
a_HOCl    = 62 * HOCl;
a_OCl     = 66 * OCl_minus;
a_NH2Cl   = 371 * NH2Cl;
a_NHCl2   = 142 * NHCl2;
a_H2O2    = 18 * H2O2;
a_HO2_rad = 228 * HO2_rad;
a_background = app.a_background;
a_sucralose      = 14.7 * sucralose;

a_total = a_HOCl + a_OCl + a_NH2Cl + a_NHCl2 + a_H2O2 + a_HO2_rad + a_background; %only for photolysis rate calculation, isn't function of time

% Quantum yields:
phi_HOCl     = 0.62;
phi_OCl_O_Cl = 0.55;
phi_OCl_O3P  = 0.074;
phi_OCl_O1D  = 0.133;
phi_NH2Cl    = 0.35;
phi_NHCl2    = 0.75;
phi_H2O2     = 0.5;
phi_HO2_rad  = 0.5;
phi_sucralose     = 0.49;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table A-5: Reactions 1 to 70
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reaction 1: NH2Cl → •NH2 + Cl•
r1 = (E0/L) * (1 - 10^(-a_total*L)) * (a_NH2Cl/a_total) * phi_NH2Cl;
dNH2Cl_dt = dNH2Cl_dt - r1;
dNH2_rad_dt = dNH2_rad_dt + r1;
dCl_rad_dt = dCl_rad_dt + r1;


% Reaction 2: NHCl2 → •NHCl + Cl•
r2 = (E0/L) * (1 - 10^(-a_total*L)) * (a_NHCl2/a_total) * phi_NHCl2;
dNHCl2_dt = dNHCl2_dt - r2;
dNHCl_rad_dt = dNHCl_rad_dt + r2;
dCl_rad_dt = dCl_rad_dt + r2;

% Reaction 3: HOCl → Cl• + •OH
r3 = (E0/L) * (1 - 10^(-a_total*L)) * (a_HOCl/a_total) * phi_HOCl;
dHOCl_dt = dHOCl_dt - r3;
dCl_rad_dt = dCl_rad_dt + r3;
dOH_rad_dt = dOH_rad_dt + r3;

% Reaction 4: NH2Cl + •OH → NHCl• + H2O
r4 = 1.02e9 * NH2Cl * OH_rad;
dNH2Cl_dt = dNH2Cl_dt - r4;
dOH_rad_dt = dOH_rad_dt - r4;
dNHCl_rad_dt = dNHCl_rad_dt + r4;

% Reaction 5: NH2Cl + Cl• → NHCl• + Cl⁻ + H⁺
r5 = 1.00e9 * NH2Cl * Cl_rad;
dNH2Cl_dt = dNH2Cl_dt - r5;
dCl_rad_dt = dCl_rad_dt - r5;
dNHCl_rad_dt = dNHCl_rad_dt + r5;
dCl_minus_dt = dCl_minus_dt + r5;
dH_plus_dt = dH_plus_dt + r5;

% Reaction 6: NH2Cl + Cl2•⁻ → NHCl• + 2Cl⁻ + H⁺
r6 = 1.14e7 * NH2Cl * Cl2_rad_minus;
dNH2Cl_dt = dNH2Cl_dt - r6;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r6;
dNHCl_rad_dt = dNHCl_rad_dt + r6;
dCl_minus_dt = dCl_minus_dt + 2*r6;
dH_plus_dt = dH_plus_dt + r6;

% Reaction 7: NHCl2 + •OH → •NCl2 + H2O
r7 = 6.21e8 * NHCl2 * OH_rad;
dNHCl2_dt = dNHCl2_dt - r7;
dOH_rad_dt = dOH_rad_dt - r7;
dNCl2_rad_dt = dNCl2_rad_dt + r7;

% Reaction 8: NHCl2 + Cl• → •NCl2 + Cl⁻ + H⁺
r8 = 1.00e9 * NHCl2 * Cl_rad;
dNHCl2_dt = dNHCl2_dt - r8;
dCl_rad_dt = dCl_rad_dt - r8;
dNCl2_rad_dt = dNCl2_rad_dt + r8;
dCl_minus_dt = dCl_minus_dt + r8;
dH_plus_dt = dH_plus_dt + r8;

% Reaction 9: NHCl2 + Cl2•⁻ → •NCl2 + 2Cl⁻ + H⁺
r9 = 4.4e6 * NHCl2 * Cl2_rad_minus;
dNHCl2_dt = dNHCl2_dt - r9;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r9;
dNCl2_rad_dt = dNCl2_rad_dt + r9;
dCl_minus_dt = dCl_minus_dt + 2*r9;
dH_plus_dt = dH_plus_dt + r9;

% Reaction 10: HOCl + •OH → ClO• + H2O
r10 = 5.00e8 * HOCl * OH_rad;
dHOCl_dt = dHOCl_dt - r10;
dOH_rad_dt = dOH_rad_dt - r10;
dClO_rad_dt = dClO_rad_dt + r10;

% Reaction 11: HOCl + Cl• → ClO• + H⁺ + Cl⁻
r11 = 3.00e9 * HOCl * Cl_rad;
dHOCl_dt = dHOCl_dt - r11;
dCl_rad_dt = dCl_rad_dt - r11;
dClO_rad_dt = dClO_rad_dt + r11;
dH_plus_dt = dH_plus_dt + r11;
dCl_minus_dt = dCl_minus_dt + r11;

% Reaction 12: OCl⁻ + •OH → ClO• + OH⁻
r12 = 1.85e9 * OCl_minus * OH_rad;
dOCl_minus_dt = dOCl_minus_dt - r12;
dOH_rad_dt = dOH_rad_dt - r12;
dClO_rad_dt = dClO_rad_dt + r12;
dOH_minus_dt = dOH_minus_dt + r12;

% Reaction 13: Cl• + OCl⁻ → ClO• + Cl⁻
r13 = 8.30e9 * Cl_rad * OCl_minus;
dCl_rad_dt = dCl_rad_dt - r13;
dOCl_minus_dt = dOCl_minus_dt - r13;
dClO_rad_dt = dClO_rad_dt + r13;
dCl_minus_dt = dCl_minus_dt + r13;

% Reaction 14: ClO• + ClO• → Cl2O2
r14 = 2.50e9 * (ClO_rad)^2;
dClO_rad_dt = dClO_rad_dt - 2*r14;
dCl2O2_dt = dCl2O2_dt + r14;

% Reaction 15: Cl2O2 + H2O → ClO2⁻ + HOCl + H⁺
r15 = 1.00e4 * Cl2O2;
dCl2O2_dt = dCl2O2_dt - r15;
dClO2_minus_dt = dClO2_minus_dt + r15;
dHOCl_dt = dHOCl_dt + r15;
dH_plus_dt = dH_plus_dt + r15;

% Reaction 16: ClO• + ClO2⁻ → ClO2· + OCl⁻
r16 = 9.40e8 * ClO_rad * ClO2_minus;
dClO_rad_dt = dClO_rad_dt - r16;
dClO2_minus_dt = dClO2_minus_dt - r16;
dClO2_rad_dt = dClO2_rad_dt + r16;
dOCl_minus_dt = dOCl_minus_dt + r16;

% Reaction 17: Cl• + ClO2⁻ → ClO2· + Cl⁻
r17 = 7.00e9 * Cl_rad * ClO2_minus;
dCl_rad_dt = dCl_rad_dt - r17;
dClO2_minus_dt = dClO2_minus_dt - r17;
dClO2_rad_dt = dClO2_rad_dt + r17;
dCl_minus_dt = dCl_minus_dt + r17;

% Reaction 18: •OH + ClO2⁻ → ClO2· + OH⁻
r18 = 7.00e9 * OH_rad * ClO2_minus;
dOH_rad_dt = dOH_rad_dt - r18;
dClO2_minus_dt = dClO2_minus_dt - r18;
dClO2_rad_dt = dClO2_rad_dt + r18;
dOH_minus_dt = dOH_minus_dt + r18;

% Reaction 19: •OH + ClO2· → ClO3⁻ + H⁺
r19 = 4.00e9 * OH_rad * ClO2_rad;
dOH_rad_dt = dOH_rad_dt - r19;
dClO2_rad_dt = dClO2_rad_dt - r19;
dClO3_minus_dt = dClO3_minus_dt + r19;
dH_plus_dt = dH_plus_dt + r19;

% Reaction 20: Cl• + ClO2· → products
r20 = 4.00e9 * Cl_rad * ClO2_rad;
dCl_rad_dt = dCl_rad_dt - r20;
dClO2_rad_dt = dClO2_rad_dt - r20;
% (Products not tracked)

% Reaction 21: •OH + ClO3⁻ → products
r21 = 1.00e6 * OH_rad * ClO3_minus;
dOH_rad_dt = dOH_rad_dt - r21;
dClO3_minus_dt = dClO3_minus_dt - r21;

% Reaction 22: Cl• + ClO3⁻ → products
r22 = 1.00e6 * Cl_rad * ClO3_minus;
dCl_rad_dt = dCl_rad_dt - r22;
dClO3_minus_dt = dClO3_minus_dt - r22;

% Reaction 23: •NH2 + O2 → NH2O2·
r23 = 1.2e8 * NH2_rad * O2;
dNH2_rad_dt = dNH2_rad_dt - r23;
dNH2O2_rad_dt = dNH2O2_rad_dt + r23;

% Reaction 24: NHCl• + O2 → NHClO2·
r24 = 1.2e8 * NHCl_rad * O2;
dNHCl_rad_dt = dNHCl_rad_dt - r24;
dNHClO2_rad_dt = dNHClO2_rad_dt + r24;

% Reaction 25: NH2O2· → •NO + H2O
r25 = 1.0e8 * NH2O2_rad;
dNH2O2_rad_dt = dNH2O2_rad_dt - r25;
dNO_rad_dt = dNO_rad_dt + r25;

% Reaction 26: NHClO2· → •NO + product
r26 = 1.0e8 * NHClO2_rad;
dNHClO2_rad_dt = dNHClO2_rad_dt - r26;
dNO_rad_dt = dNO_rad_dt + r26;

% Reaction 27: NH2O2· → transient species → N2O
r27 = 5.98e8 * NH2O2_rad;
dNH2O2_rad_dt = dNH2O2_rad_dt - r27;
dN2O_dt = dN2O_dt + r27;

% Reaction 28: NHClO2· → transient species → N2O
r28 = 6.70e8 * NHClO2_rad;
dNHClO2_rad_dt = dNHClO2_rad_dt - r28;
dN2O_dt = dN2O_dt + r28;

% Reaction 29: •NO + •OH → NO2⁻ + H⁺
r29 = 1.0e10 * NO_rad * OH_rad;
dNO_rad_dt = dNO_rad_dt - r29;
dOH_rad_dt = dOH_rad_dt - r29;
dNO2_minus_dt = dNO2_minus_dt + r29;
dH_plus_dt = dH_plus_dt + r29;

% Reaction 30: NO2⁻ + •OH → •NO2 + OH⁻
r30 = 1.2e10 * NO2_minus * OH_rad;
dNO2_minus_dt = dNO2_minus_dt - r30;
dOH_rad_dt = dOH_rad_dt - r30;
dNO2_rad_dt = dNO2_rad_dt + r30;
dOH_minus_dt = dOH_minus_dt + r30;

% Reaction 31: •NO + •NO + O2 → 2•NO2
r31 = 2.1e6 * (NO_rad)^2 * O2;
dNO_rad_dt = dNO_rad_dt - 2*r31;
dNO2_rad_dt = dNO2_rad_dt + 2*r31;

% Reaction 32: NO + •NO2 → N2O3
r32 = 1.1e9 * NO_rad * NO2_rad;
dNO_rad_dt = dNO_rad_dt - r32;
dNO2_rad_dt = dNO2_rad_dt - r32;
dN2O3_dt = dN2O3_dt + r32;

% Reaction 33: N2O3 → •NO + •NO2
r33 = 4.3e6 * N2O3;
dN2O3_dt = dN2O3_dt - r33;
dNO_rad_dt = dNO_rad_dt + r33;
dNO2_rad_dt = dNO2_rad_dt + r33;

% Reaction 34: N2O3 + H2O → 2NO2⁻ + 2H⁺
r34 = 1.6e3 * N2O3;
dN2O3_dt = dN2O3_dt - r34;
dNO2_minus_dt = dNO2_minus_dt + 2*r34;
dH_plus_dt = dH_plus_dt + 2*r34;

% Reaction 35: •NO2 + •NO2 → N2O4
r35 = 4.5e8 * (NO2_rad)^2;
dNO2_rad_dt = dNO2_rad_dt - 2*r35;
dN2O4_dt = dN2O4_dt + r35;

% Reaction 36: N2O4 + H2O → NO2⁻ + NO3⁻ + 2H⁺
r36 = 1.0e3 * N2O4;
dN2O4_dt = dN2O4_dt - r36;
dNO2_minus_dt = dNO2_minus_dt + r36;
dNO3_minus_dt = dNO3_minus_dt + r36;
dH_plus_dt = dH_plus_dt + 2*r36;

% Reaction 37 (Table A-5, Reaction 37): H2O → H⁺ + OH⁻
r37_extra = 1.00e-3 ; % assumed 55.5 moles/L of H2O (maybe it's built in?)
dH_plus_dt = dH_plus_dt + r37_extra;
dOH_minus_dt = dOH_minus_dt + r37_extra;

% Reaction 38 (Table A-5, Reaction 38): H⁺ + OH⁻ → H2O
r38_extra = 1.00e11 * H_plus * OH_minus;
dH_plus_dt = dH_plus_dt - r38_extra;
dOH_minus_dt = dOH_minus_dt - r38_extra;

% Reaction 39 (Table A-5, Reaction 39): •OH + •OH → H2O2
r39_extra = 5.50e9 * (OH_rad)^2;
dOH_rad_dt = dOH_rad_dt - 2*r39_extra;
dH2O2_dt = dH2O2_dt + r39_extra;

% Reaction 40 (Table A-5, Reaction 40): •OH + OH⁻ → O•⁻ + H2O
r40_extra = 1.20e10 * OH_rad * OH_minus;
dOH_rad_dt = dOH_rad_dt - r40_extra;
dO_rad_minus_dt = dO_rad_minus_dt + r40_extra;  % Now tracked as species 49

% Reaction 41 (Table A-5, Reaction 41): •OH + Cl⁻ → ClOH•⁻
r41_extra = 4.30e9 * OH_rad * Cl_minus;
dOH_rad_dt = dOH_rad_dt - r41_extra;
dClOH_rad_minus_dt = dClOH_rad_minus_dt + r41_extra;
dCl_minus_dt = dCl_minus_dt - r41_extra;

% Reaction 42 (Table A-5, Reaction 42): ClOH•⁻ → Cl⁻ + •OH
r42_extra = 6.10e9 * ClOH_rad_minus;
dClOH_rad_minus_dt = dClOH_rad_minus_dt - r42_extra;
dCl_minus_dt = dCl_minus_dt + r42_extra;
dOH_rad_dt = dOH_rad_dt + r42_extra;

% Reaction 43 (Table A-5, Reaction 43): Cl• + H2O → ClOH•⁻ + H⁺
r43_extra = 2.50e5 * Cl_rad * 55.5;
dCl_rad_dt = dCl_rad_dt - r43_extra;
dClOH_rad_minus_dt = dClOH_rad_minus_dt + r43_extra;
dH_plus_dt = dH_plus_dt + r43_extra;

% Reaction 44 (Table A-5, Reaction 44): Cl• + OH⁻ → ClOH•⁻
r44_extra = 1.80e10 * Cl_rad * OH_minus;
dCl_rad_dt = dCl_rad_dt - r44_extra;
dOH_minus_dt = dOH_minus_dt - r44_extra;
dClOH_rad_minus_dt = dClOH_rad_minus_dt + r44_extra;

% Reaction 45 (Table A-5, Reaction 45): Cl• + Cl⁻ → Cl2•⁻
r45_extra = 8.00e9 * Cl_rad * Cl_minus;
dCl_rad_dt = dCl_rad_dt - r45_extra;
dCl_minus_dt = dCl_minus_dt - r45_extra;
dCl2_rad_minus_dt = dCl2_rad_minus_dt + r45_extra;

% Reaction 46 (Table A-5, Reaction 46): Cl2•⁻ → Cl• + Cl⁻
r46_extra = 6.00e4 * Cl2_rad_minus;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r46_extra;
dCl_rad_dt = dCl_rad_dt + r46_extra;
dCl_minus_dt = dCl_minus_dt + r46_extra;

% Reaction 47 (Table A-5, Reaction 47): Cl2•⁻ + H2O → Cl⁻ + HClOH•
r47_extra = 1.30e3 * Cl2_rad_minus * 55.5;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r47_extra;
dCl_minus_dt = dCl_minus_dt + r47_extra;
dHClOH_rad_dt = dHClOH_rad_dt + r47_extra;

% Reaction 48 (Table A-5, Reaction 48): HClOH• → ClOH•⁻ + H⁺
r48_extra = 1.00e8 * HClOH_rad;
dHClOH_rad_dt = dHClOH_rad_dt - r48_extra;
dClOH_rad_minus_dt = dClOH_rad_minus_dt + r48_extra;
dH_plus_dt = dH_plus_dt + r48_extra;

% Reaction 49 (Table A-5, Reaction 49): ClOH•⁻ + Cl⁻ → Cl2•⁻ + OH⁻
r49_extra = 1.00e4 * ClOH_rad_minus * Cl_minus;
dClOH_rad_minus_dt = dClOH_rad_minus_dt - r49_extra;
dCl_minus_dt = dCl_minus_dt - r49_extra;
dCl2_rad_minus_dt = dCl2_rad_minus_dt + r49_extra;
dOH_minus_dt = dOH_minus_dt + 2*r49_extra;

% Reaction 50 (Table A-5, Reaction 50): ClOH•⁻ + H⁺ → Cl• + H2O
r50_extra = 2.10e10 * ClOH_rad_minus * H_plus;
dClOH_rad_minus_dt = dClOH_rad_minus_dt - r50_extra;
dH_plus_dt = dH_plus_dt - r50_extra;
dCl_rad_dt = dCl_rad_dt + r50_extra;

% Reaction 51 (Table A-5, Reaction 51): Cl• + Cl• → Cl2
r51_extra = 8.80e7 * (Cl_rad)^2;
dCl_rad_dt = dCl_rad_dt - 2*r51_extra;
dCl2_dt = dCl2_dt + r51_extra;

% Reaction 52 (Table A-5, Reaction 52): Cl2•⁻ + •OH → HOCl + Cl⁻
r52_extra = 1.00e9 * Cl2_rad_minus * OH_rad;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r52_extra;
dOH_rad_dt = dOH_rad_dt - r52_extra;
dHOCl_dt = dHOCl_dt + r52_extra;
dCl_minus_dt = dCl_minus_dt + r52_extra;

% Reaction 53 (Table A-5, Reaction 53): Cl2•⁻ + Cl2•⁻ → Cl2 + 2Cl⁻
r53_extra = 9.00e8 * (Cl2_rad_minus)^2;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r53_extra;
dCl2_dt = dCl2_dt + r53_extra;
dCl_minus_dt = dCl_minus_dt + 2*r53_extra;

% Reaction 54 (Table A-5, Reaction 54): Cl2•⁻ + Cl• → Cl2 + Cl⁻
r54_extra = 2.10e9 * Cl2_rad_minus * Cl_rad;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r54_extra;
dCl_rad_dt = dCl_rad_dt - r54_extra;
dCl2_dt = dCl2_dt + r54_extra;
dCl_minus_dt = dCl_minus_dt + r54_extra;

% Reaction 55 (Table A-5, Reaction 55): Cl2•⁻ + OH⁻ → ClOH•⁻ + Cl⁻
r55_extra = 4.50e7 * Cl2_rad_minus * OH_minus;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r55_extra;
dOH_minus_dt = dOH_minus_dt - r55_extra;
dClOH_rad_minus_dt = dClOH_rad_minus_dt + r55_extra;
dCl_minus_dt = dCl_minus_dt + r55_extra;

% Reaction 56 (Table A-5, Reaction 56): Cl3⁻ → Cl2 + Cl⁻
r56_extra = 1.10e5 * Cl3_minus;
dCl3_minus_dt = dCl3_minus_dt - r56_extra;
dCl2_dt = dCl2_dt + r56_extra;
dCl_minus_dt = dCl_minus_dt + r56_extra;

% Reaction 57 (Table A-5, Reaction 57): Cl2 + H2O → HOCl + Cl⁻ + H⁺
r57_extra = 15 * Cl2;
dCl2_dt = dCl2_dt - r57_extra;
dHOCl_dt = dHOCl_dt + r57_extra;
dCl_minus_dt = dCl_minus_dt + r57_extra;
dH_plus_dt = dH_plus_dt + r57_extra;

% Reaction 58 (Table A-5, Reaction 58): Cl2 + Cl⁻ → Cl3⁻
r58_extra = 2.00e4 * Cl2 * Cl_minus;
dCl2_dt = dCl2_dt - r58_extra;
dCl_minus_dt = dCl_minus_dt - r58_extra;
dCl3_minus_dt = dCl3_minus_dt + r58_extra;

% Reaction 59 (Table A-5, Reaction 59): HOCl → OCl⁻ + H⁺
r59_extra = 1.41e3 * HOCl;
dHOCl_dt = dHOCl_dt - r59_extra;
dOCl_minus_dt = dOCl_minus_dt + r59_extra;
dH_plus_dt = dH_plus_dt + r59_extra;

% Reaction 60 (Table A-5, Reaction 60): OCl⁻ + H⁺ → HOCl
r60_extra = 5.00e10 * OCl_minus * H_plus;
dOCl_minus_dt = dOCl_minus_dt - r60_extra;
dH_plus_dt = dH_plus_dt - r60_extra;
dHOCl_dt = dHOCl_dt + r60_extra;

% Reaction 61 (Table A-5, Reaction 61): H⁺ + Cl⁻ → HCl
r61_extra = 5.00e10 * H_plus * Cl_minus;
dH_plus_dt = dH_plus_dt - r61_extra;
dCl_minus_dt = dCl_minus_dt - r61_extra;
dHCl_dt = dHCl_dt + r61_extra;

% Reaction 62 (Table A-5, Reaction 62): HCl → H⁺ + Cl⁻
r62_extra = 8.60e16 * HCl;
dHCl_dt = dHCl_dt - r62_extra;
dH_plus_dt = dH_plus_dt + r62_extra;
dCl_minus_dt = dCl_minus_dt + r62_extra;

% Reaction 63 (Table A-5, Reaction 63): H⁺ + NH2Cl + NO2⁻ → NH3 + NO2Cl → NO3⁻
r63_extra = 1.36e7 * H_plus * NH2Cl * NO2_minus;
dH_plus_dt = dH_plus_dt - r63_extra;
dNH2Cl_dt = dNH2Cl_dt - r63_extra;
dNO2_minus_dt = dNO2_minus_dt - r63_extra;
dNH3_dt = dNH3_dt + r63_extra;
dNO3_minus_dt = dNO3_minus_dt + r63_extra;

% Reaction 64 (Table A-5, Reaction 64): •OH + HCO3⁻ → CO3·⁻ + H2O
r64_extra = 8.50e6 * OH_rad * HCO3_minus;
dOH_rad_dt = dOH_rad_dt - r64_extra;
dHCO3_minus_dt = dHCO3_minus_dt - r64_extra;
dCO3_rad_dt = dCO3_rad_dt + r64_extra;

% Reaction 65 (Table A-5, Reaction 65): •OH + H2CO3 → CO3·⁻ + H2O + H⁺
r65_extra = 1.00e6 * OH_rad * H2CO3;
dOH_rad_dt = dOH_rad_dt - r65_extra;
dH2CO3_dt = dH2CO3_dt - r65_extra;
dCO3_rad_dt = dCO3_rad_dt + r65_extra;
dH_plus_dt = dH_plus_dt + r65_extra;

% Reaction 66 (Table A-5, Reaction 66): •OH + CO3·⁻ → products
r66_extra = 3.00e9 * OH_rad * CO3_rad;
dOH_rad_dt = dOH_rad_dt - r66_extra;
dCO3_rad_dt = dCO3_rad_dt - r66_extra;

% Reaction 67 (Table A-5, Reaction 67): CO3·⁻ + CO3·⁻ → products
r67_extra = 3.00e7 * (CO3_rad)^2;
dCO3_rad_dt = dCO3_rad_dt - r67_extra;

% Reaction 68 (Table A-5, Reaction 68): 1,4-dioxane + •OH → products
r68_extra = 3.1e9 * dioxane * OH_rad;
ddioxane_dt = ddioxane_dt - r68_extra;
dOH_rad_dt = dOH_rad_dt - r68_extra;

% Reaction 69 (Table A-5, Reaction 69): 1,4-dioxane + Cl• → products
r69_extra = 4.4e6 * dioxane * Cl_rad;
ddioxane_dt = ddioxane_dt - r69_extra;
dCl_rad_dt = dCl_rad_dt - r69_extra;

% Reaction 70 (Table A-5, Reaction 70): 1,4-dioxane + Cl2•⁻ → products
r70_extra = 3.3e6 * dioxane * Cl2_rad_minus;
ddioxane_dt = ddioxane_dt - r70_extra;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r70_extra;

% Reaction (report, Reaction caffeine_OH): caffeine + •OH → products
rcaffeine_OH = 8.8e9 * caffeine * OH_rad;
dcaffeine_dt = dcaffeine_dt - rcaffeine_OH;
dOH_rad_dt = dOH_rad_dt - rcaffeine_OH;

% Reaction (report, Reaction caffeine_Cl): caffeine + Cl• → products
rcaffeine_Cl = 3.3e10 * caffeine * Cl_rad;
dcaffeine_dt = dcaffeine_dt - rcaffeine_Cl;
dCl_rad_dt = dCl_rad_dt - rcaffeine_Cl;
% 
% Reaction (report, Reaction caffeine_Cl2): caffeine + Cl2•⁻ → products
rcaffeine_Cl2 = 0.43e8 * caffeine * Cl2_rad_minus;
dcaffeine_dt = dcaffeine_dt - rcaffeine_Cl2;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - rcaffeine_Cl2;

% Reaction (report, Reaction caffeine_ClO): caffeine + ClO•⁻ → products
rcaffeine_Cl2 = 1.03e8 * caffeine * ClO_rad;
dcaffeine_dt = dcaffeine_dt - rcaffeine_Cl2;
dClO_rad_dt = dClO_rad_dt  - rcaffeine_Cl2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reaction (report, Reaction sucralose_OH): sucralose + •OH → products
rsucralose_OH = 1.6e9 * sucralose * OH_rad;
dsucralose_dt = dsucralose_dt - rsucralose_OH;
dOH_rad_dt = dOH_rad_dt - rcaffeine_OH;

% Reaction (report, Reaction sucralose_Cl): sucralose+ Cl• → products
rsucralose_Cl = 1.11e10 * sucralose * Cl_rad;
dsucralose_dt = dsucralose_dt - rsucralose_Cl;
dCl_rad_dt = dCl_rad_dt - rsucralose_Cl;
% 
% Reaction (report, Reaction sucralose_Cl2): sucralose+ Cl2•⁻ → products
rsucralose_Cl2 = 0.01e8 * sucralose * Cl2_rad_minus;
dsucralose_dt = dsucralose_dt - rsucralose_Cl2;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - rsucralose_Cl2;


%%%%%%%%%%%%%%%%%%%%
% %% Table A-1 Reactions (Model numbering from the table)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Reaction 1: HOCl + NH3 → NH2Cl + H2O
r71 = 2.04e9 * exp(-1887/T) * HOCl * NH3;
dHOCl_dt = dHOCl_dt - r71;
dNH3_dt  = dNH3_dt - r71;
dNH2Cl_dt = dNH2Cl_dt + r71;

% Reaction 2: NH2Cl + H2O → HOCl + NH3
r72 = 1.38e8 * exp(-8800/T) * NH2Cl;
dNH2Cl_dt = dNH2Cl_dt - r72;
dHOCl_dt = dHOCl_dt + r72;
dNH3_dt  = dNH3_dt + r72;

% Reaction 3: HOCl + NH2Cl → NHCl2 + H2O
r73 = 3.0e5 * exp(-2010/T) * HOCl * NH2Cl;
dHOCl_dt = dHOCl_dt - r73;
dNH2Cl_dt = dNH2Cl_dt - r73;
dNHCl2_dt = dNHCl2_dt + r73;

% Reaction 4: NHCl2 + H2O → HOCl + NH2Cl
r74 = 6.5e-7 * NHCl2;
dNHCl2_dt = dNHCl2_dt - r74;
dHOCl_dt = dHOCl_dt + r74;
dNH2Cl_dt = dNH2Cl_dt + r74;

% Reaction 5: NH2Cl + NH2Cl → NHCl2 + NH3
kH    = 3.78e10 * exp(-2169/T) / 3600;
kHCO3 = 0.87 * exp(-503/T) / 3600;
kH2CO3 = 2.52e25 * exp(-16860/T) / 3600;
r75 = (kH * H_plus + kHCO3 * HCO3_minus + kH2CO3 * H2CO3) * (NH2Cl)^2;
dNH2Cl_dt = dNH2Cl_dt - r75;
dNHCl2_dt = dNHCl2_dt + r75;
dNH3_dt  = dNH3_dt + r75;

% Reaction 6: NHCl2 + NH3 → 2 NH2Cl
r76 = 2.67e4 * NHCl2 * NH3;
dNHCl2_dt = dNHCl2_dt - r76;
dNH3_dt   = dNH3_dt - r76;
dNH2Cl_dt = dNH2Cl_dt + 2*r76;

% Reaction 7: NHCl2 + H2O → NOH + 2H⁺ + 2Cl⁻
r77 = 1.67e2 * NHCl2 * OH_minus;  % using [OH⁻] = OH_minus
dNHCl2_dt = dNHCl2_dt - r77;
dNOH_dt   = dNOH_dt + r77;
dH_plus_dt = dH_plus_dt + 2*r77;
dCl_minus_dt = dCl_minus_dt + 2*r77;

% Reaction 8: NOH + NHCl2 → HOCl + N2 + H⁺ + Cl⁻
r78 = 2.77e4 * NOH * NHCl2;
dNOH_dt   = dNOH_dt - r78;
dNHCl2_dt = dNHCl2_dt - r78;
dHOCl_dt  = dHOCl_dt + r78;
dN2_dt    = dN2_dt + r78;
dH_plus_dt = dH_plus_dt + r78;
dCl_minus_dt = dCl_minus_dt + r78;

% Reaction 9: NOH + NH2Cl → N2 + HCl + H2O
r79 = 8.3e3 * NOH * NH2Cl;
dNOH_dt  = dNOH_dt - r79;
dNH2Cl_dt = dNH2Cl_dt - r79;
dN2_dt   = dN2_dt + r79;
dHCl_dt  = dHCl_dt + r79;

% Reaction 10: NH2Cl + NHCl2 → N2 + 3H + 3Cl⁻ (k10 = 0)
r80 = 0;

% Reaction 11: HOCl + NHCl2 → NCl3 + H2O
k11 = 3.28e9 * OH_minus + 9.00e4 * OCl_minus + 6.00e6 * CO3_2;  % using [OH⁻]=OH_minus, [OCl⁻]=OCl_minus, [CO3²⁻]=CO3_2
r81 = k11 * HOCl * NHCl2;
dHOCl_dt = dHOCl_dt - r81;
dNHCl2_dt = dNHCl2_dt - r81;
dNCl3_dt = dNCl3_dt + r81;

% Reaction 12: NHCl2 + NCl3 + 2H2O → 2HOCl + N2 + 3H⁺ + 3Cl⁻
r82 = 1.00e14 * NHCl2 * NCl3;
dNHCl2_dt = dNHCl2_dt - r82;
dNCl3_dt = dNCl3_dt - r82;
dHOCl_dt = dHOCl_dt + 2*r82;
dN2_dt   = dN2_dt + r82;
dH_plus_dt = dH_plus_dt + 3*r82;
dCl_minus_dt = dCl_minus_dt + 3*r82;

% Reaction 13: NH2Cl + NCl3 + H2O → HOCl + N2 + 3H⁺ + 3Cl⁻
r83 = 1.00e6 * NH2Cl * NCl3;
dNH2Cl_dt = dNH2Cl_dt - r83;
dNCl3_dt = dNCl3_dt - r83;
dHOCl_dt = dHOCl_dt + r83;
dN2_dt   = dN2_dt + r83;
dH_plus_dt = dH_plus_dt + 3*r83;
dCl_minus_dt = dCl_minus_dt + 3*r83;

% Reaction 14: NHCl2 + 2HOCl + H2O → NO3⁻ + 5H⁺ + 4Cl⁻
r84 = 66 * NHCl2 * OCl_minus;  % using OCl_minus as proxy for [OCl⁻]
dNHCl2_dt = dNHCl2_dt - r84;
dHOCl_dt  = dHOCl_dt - r84;
dNO3_minus_dt = dNO3_minus_dt + r84;
dH_plus_dt = dH_plus_dt + 5*r84;
dCl_minus_dt = dCl_minus_dt + 4*r84;

% Reaction 15: NCl3 + H2O → NHCl2 + HOCl
r85 = (1.60e-6 + 8*OH_minus + 890*(OH_minus)^2 + 65*OH_minus*HCO3_minus) * NCl3;
dNCl3_dt = dNCl3_dt - r85;
dNHCl2_dt = dNHCl2_dt + r85;
dHOCl_dt = dHOCl_dt + r85;

% Reaction 16: (Not defined)
r86 = 0;

% Reaction 17: (Not defined)
r87 = 0;

% Reaction 18: NH4⁺ → NH3 + H⁺
r88 = 42 * NH4;
dNH4_dt = dNH4_dt - r88;
dNH3_dt = dNH3_dt + r88;
dH_plus_dt = dH_plus_dt + r88;

% Reaction 19: NH3 + H⁺ → NH4⁺
r89 = 8.4e10 * NH3 * H_plus;
dNH3_dt = dNH3_dt - r89;
dH_plus_dt = dH_plus_dt - r89;
dNH4_dt = dNH4_dt + r89;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =======================================================
%% TABLE C-1: Oxidant Equilibria, Initiation, and Radical Reactions
%% =======================================================

% Oxidant Equilibria: 
% Reaction 1 (C-1): Cl₂ + H₂O ↔ HOCl + Cl⁻ + H⁺ (duplicate)
% 
% 
% Reaction 2 (C-1): H₂O₂ ↔ H⁺ + HO2⁻ (duplicate)
% 
% 
% Reaction 3 (C-1): HCO3⁻ ↔ CO3²⁻ + H⁺ (duplicated in 164-167 )


% --- Initiation (Photolysis) ---
% Reaction 5 (C-1): H₂O₂ + hv → 2·OH
r_C5 = (E0/L) * (1 - 10^(-a_total*L)) * (a_H2O2/a_total) * phi_H2O2;
dH2O2_dt = dH2O2_dt - r_C5;
dOH_rad_dt = dOH_rad_dt + 2*r_C5;

% Reaction 6 (C-1): HO2⁻ + H₂O + hv → 2·OH + OH⁻
r_C6 = (E0/L) * (1 - 10^(-a_total*L)) * (a_HO2_rad/a_total) * phi_HO2_rad;
dHO2_minus_dt = dHO2_minus_dt - r_C6;
dOH_rad_dt = dOH_rad_dt + 2*r_C6;
dOH_minus_dt = dOH_minus_dt + r_C6;

% Reaction 7 (C-1): HOCl + hv → ·OH + Cl• (duplicate)


% Reaction 8 (C-1): OCl⁻ + H₂O + hv → Cl• + ·OH + OH⁻ (duplicate)


% --- Radical Reactions ---
% Reaction 11 (C-1): ·OH + DOC → products, k = 3.0e8 M⁻¹ s⁻¹
r_C11 = 3.0e8 * OH_rad * DOC;
dOH_rad_dt = dOH_rad_dt - r_C11;
dDOC_dt = dDOC_dt - r_C11;

% Reaction 13 (C-1): ·OH + CO3²⁻ → CO3·⁻ + OH⁻, k = 3.9e8 M⁻¹ s⁻¹
r_C13 = 3.9e8 * OH_rad * CO3_2;
dOH_rad_dt = dOH_rad_dt - r_C13;
dCO3_rad_dt = dCO3_rad_dt + r_C13;
dOH_minus_dt = dOH_minus_dt + r_C13;

% Reaction 16 (C-1): ·OH + O•⁻ → HO2⁻, k = 1.0e10 M⁻¹ s⁻¹
r_C16 = 1.0e10 * OH_rad * O_rad_minus;
dOH_rad_dt = dOH_rad_dt - r_C16;
dO_rad_minus_dt = dO_rad_minus_dt - r_C16;
dHO2_minus_dt = dHO2_minus_dt + r_C16;

% Reaction 17 (C-1): ·OH + HO2⁻ → HO2· + OH⁻, k = 7.5e9 M⁻¹ s⁻¹
r_C17 = 7.5e9 * OH_rad * HO2_minus;
dOH_rad_dt = dOH_rad_dt - r_C17;
dHO2_rad_dt = dHO2_rad_dt + r_C17;
dOH_minus_dt = dOH_minus_dt + r_C17;

% Reaction 18 (C-1): ·OH + HO2· → H2O + O2, k = 6.6e9 M⁻¹ s⁻¹ (duplicate)


% Reaction 19 (C-1): ·OH + H2O2 → HO2· + H2O, k = 3.2e7 M⁻¹ s⁻¹
r_C19 = 3.2e7 * OH_rad * H2O2;
dOH_rad_dt = dOH_rad_dt - r_C19;
dH2O2_dt = dH2O2_dt - r_C19;
dHO2_rad_dt = dHO2_rad_dt + r_C19;

% Reaction 21 (C-1): ·OH + ClO2⁻ → ClO2· + OH⁻, k = 6.3e9 M⁻¹ s⁻¹ (duplicate)


% Reaction 23 (C-1): ·OH + O₂·⁻ → OH⁻ + O₂, k = 7.0e9 M⁻¹ s⁻¹
r_C23 = 7.0e9 * OH_rad * O2_super;
dOH_rad_dt = dOH_rad_dt - r_C23;
dO2_super_dt = dO2_super_dt - r_C23;
dOH_minus_dt = dOH_minus_dt + r_C23;
dO2_dt = dO2_dt + r_C23;

% Reaction  24 (C-1): •OH + HPO₄²⁻ → HPO₄•⁻ + OH⁻ | replaced by k= 1.5e5 (Chuang et al 2017)
% r_C1_24 = 1.5e9 * OH_rad * HPO4_2;
r_C1_24 = 1.5e5 * OH_rad * HPO4_2;
dOH_rad_dt   = dOH_rad_dt   - r_C1_24;   % •OH consumed
dHPO4_2_dt   = dHPO4_2_dt   - r_C1_24;   % HPO4^2- consumed
dHPO4_rad_dt = dHPO4_rad_dt + r_C1_24;   % HPO4•- produced
dOH_minus_dt = dOH_minus_dt + r_C1_24;     % OH^- produced

% Reaction 25 (C-1): •OH + H2PO₄1⁻ → HPO₄•⁻ + OH⁻ | replaced by k= 2e4 (Chuang et al 2017)
% r_C1_25 = 2.0e9 * OH_rad * H2PO4_minus;
r_C1_25 = 2.0e4 * OH_rad * H2PO4_minus;
dOH_rad_dt      = dOH_rad_dt      - r_C1_25;  % •OH consumed
dH2PO4_minus_dt = dH2PO4_minus_dt - r_C1_25;  % H2PO4^- consumed
dHPO4_rad_dt    = dHPO4_rad_dt    + r_C1_25;  % HPO4•- produced


% Reaction 28 (C-1): Cl• + HCO3⁻ → CO3·⁻ + Cl⁻ + H⁺, k = 2.2e8 M⁻¹ s⁻¹
r_C28 = 2.2e8 * Cl_rad * HCO3_minus;
dCl_rad_dt = dCl_rad_dt - r_C28;
dCO3_rad_dt = dCO3_rad_dt + r_C28;
dCl_minus_dt = dCl_minus_dt + r_C28;
dH_plus_dt = dH_plus_dt + r_C28;

% Reaction 29 (C-1): Cl• + CO3²⁻ → CO3·⁻ + Cl⁻, k = 5.0e8 M⁻¹ s⁻¹
r_C29 = 5.0e8 * Cl_rad * CO3_2;
dCl_rad_dt = dCl_rad_dt - r_C29;
dCO3_rad_dt = dCO3_rad_dt + r_C29;
dCl_minus_dt = dCl_minus_dt + r_C29;

% Reaction 32 (C-1): Cl• + Cl2•⁻ → Cl⁻ + Cl₂, k = 2.1e9 M⁻¹ s⁻¹
r_C32 = 2.1e9 * Cl_rad * Cl2_rad_minus;
dCl_rad_dt = dCl_rad_dt - r_C32;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r_C32;
dCl_minus_dt = dCl_minus_dt + r_C32;
dCl2_dt = dCl2_dt + r_C32;

% Reaction 35 (C-1): Cl• + DOC → products, k = 1.56e8 M⁻¹ s⁻¹
r_C35 = 1.56e8 * Cl_rad * DOC;
dCl_rad_dt = dCl_rad_dt - r_C35;
dDOC_dt = dDOC_dt - r_C35;

% Reaction 36 (C-1): Cl• + H2O2 → H⁺ + Cl⁻ + HO2•, k = 2.0e9 M⁻¹ s⁻¹
r_C36 = 2.0e9 * Cl_rad * H2O2;
dCl_rad_dt = dCl_rad_dt - r_C36;
dH2O2_dt = dH2O2_dt - r_C36;
dH_plus_dt = dH_plus_dt + r_C36;
dCl_minus_dt = dCl_minus_dt + r_C36;
dHO2_rad_dt = dHO2_rad_dt + r_C36;

% --- Reaction 41 ---
% Cl₂·⁻ + HCO₃⁻ → CO₃·⁻ + 2Cl⁻ + H⁺      k = 8.0×10^7 M⁻¹·s⁻¹
r41 = 8.0e7 * Cl2_rad_minus * HCO3_minus;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r41;
dCO3_rad_dt       = dCO3_rad_dt       + r41;
dCl_minus_dt      = dCl_minus_dt      + 2*r41;
dH_plus_dt        = dH_plus_dt        + r41;

% --- Reaction 42 ---
% (No reaction provided; nothing to code)

% --- Reaction 43 ---
% (No reaction provided)

% --- Reaction 44 ---
% (No reaction provided)

% --- Reaction 45 ---
% Cl₂·⁻ + H₂O → Cl⁻ + HClOH·      k = 1.3×10^3 s⁻¹
% (Here the water term is assumed constant and absorbed into the rate constant.)
r45 = 1.3e3 * Cl2_rad_minus;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r45;
dCl_minus_dt      = dCl_minus_dt      + r45;
dHClOH_rad_dt     = dHClOH_rad_dt     + r45;

% --- Reaction 46 ---
% (No reaction provided)

% --- Reaction 47 ---
% (No reaction provided)

% --- Reaction 48 ---
% (No reaction provided)

% --- Reaction 49 ---
% Cl₂·⁻ + ClO₂⁻ → ClO₂· + 2Cl⁻      k = 1.3×10^8 M⁻¹·s⁻¹
r49 = 1.3e8 * Cl2_rad_minus * ClO2_minus;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r49;
dClO2_minus_dt    = dClO2_minus_dt    - r49;
dClO2_rad_dt      = dClO2_rad_dt      + r49;
dCl_minus_dt      = dCl_minus_dt      + 2*r49;

% --- Reaction 50 ---
% Cl₂·⁻ + CO₃²⁻ → CO₃·⁻ + 2Cl⁻      k = 1.6×10^8 M⁻¹·s⁻¹
r50 = 1.6e8 * Cl2_rad_minus * CO3_2;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r50;
dCO3_rad_dt       = dCO3_rad_dt       + r50;
dCl_minus_dt      = dCl_minus_dt      + 2*r50;

% --- Reaction 51 ---
% (No reaction provided)

% --- Reaction 52 ---
% ClO· + DOC → products      k = 1.56×10^8 M⁻¹·s⁻¹
r52 = 1.56e8 * ClO_rad * DOC;
dClO_rad_dt = dClO_rad_dt - r52;
dDOC_dt     = dDOC_dt     - r52;
% (Products are not tracked)

% --- Reaction 53 ---
% (No reaction provided)

% --- Reaction 54 ---
% (No reaction provided)

% --- Reaction 55 ---
% 2ClO· + H₂O → HOCl + H⁺ + ClO₂⁻      k = 2.5×10^9 M⁻¹·s⁻¹
% (Assuming water is in excess so that its concentration is absorbed into the rate constant)
r55 = 2.5e9 * (ClO_rad)^2;
dClO_rad_dt    = dClO_rad_dt    - 2*r55;
dHOCl_dt       = dHOCl_dt       + r55;
dH_plus_dt     = dH_plus_dt     + r55;
dClO2_minus_dt = dClO2_minus_dt + r55;

% --- Reaction 56 ---
% (No reaction provided)

% --- Reaction 57 ---
% CO₃·⁻ + ·OH → products      k = 3.0×10^9 M⁻¹·s⁻¹
r57 = 3.0e9 * CO3_rad * OH_rad;
dCO3_rad_dt = dCO3_rad_dt - r57;
dOH_rad_dt  = dOH_rad_dt  - r57;
% (Products not tracked)

% --- Reaction 58 ---
% (No reaction provided)

% --- Reaction 59 ---
% CO₃·⁻ + CO₃·⁻ → products      k = 3.0×10^7 M⁻¹·s⁻¹
r59 = 3.0e7 * (CO3_rad)^2;
dCO3_rad_dt = dCO3_rad_dt - 2*r59;
% (Products not tracked)

% --- Reaction 60 ---
% (No reaction provided)

% --- Reaction 61 ---
% (No reaction provided)

% --- Reaction 62 ---
% CO₃·⁻ + HO₂⁻ → HCO₃⁻ + O₂·⁻      k = 3.0×10^7 M⁻¹·s⁻¹
r62 = 3.0e7 * CO3_rad * HO2_minus;
dCO3_rad_dt    = dCO3_rad_dt    - r62;
dHO2_minus_dt  = dHO2_minus_dt  - r62;
dHCO3_minus_dt = dHCO3_minus_dt + r62;
dO2_super_dt   = dO2_super_dt   + r62;

% --- Reaction 63 ---
% (No reaction provided)

% --- Reaction 64 ---
% CO₃·⁻ + ClO₂⁻ → CO₃²⁻ + ClO₂      k = 3.4×10^7 M⁻¹·s⁻¹
r64 = 3.4e7 * CO3_rad * ClO2_minus;
dCO3_rad_dt     = dCO3_rad_dt     - r64;
dClO2_minus_dt  = dClO2_minus_dt  - r64;
dCO3_2_dt       = dCO3_2_dt       + r64;
dClO2_dt        = dClO2_dt        + r64;

% --- Reaction 65 ---
% (No reaction provided)

% --- Reaction 66 ---
% HClOH· → Cl· + H₂O      k = 1.0×10^2 s⁻¹
r66 = 1.0e2 * HClOH_rad;
dHClOH_rad_dt = dHClOH_rad_dt - r66;
dCl_rad_dt    = dCl_rad_dt    + r66;
% (H₂O is not tracked)

% --- Reaction 67 ---
% HClOH· + Cl⁻ → Cl₂·⁻ + H₂O      k = 5.0×10^9 M⁻¹·s⁻¹
r67 = 5.0e9 * HClOH_rad * Cl_minus;
dHClOH_rad_dt   = dHClOH_rad_dt   - r67;
dCl_minus_dt    = dCl_minus_dt    - r67;
dCl2_rad_minus_dt = dCl2_rad_minus_dt + r67;

% --- Reaction 68 ---
% O·⁻ + HO₂⁻ → O₂·⁻ + OH⁻      k ≈ 1.0×10^8 M⁻¹·s⁻¹ (value assumed)
k68 = 1.0e8;
r68 = k68 * O_rad_minus * HO2_minus;
dO_rad_minus_dt = dO_rad_minus_dt - r68;
dHO2_minus_dt   = dHO2_minus_dt   - r68;
dO2_super_dt    = dO2_super_dt    + r68;
dOH_minus_dt    = dOH_minus_dt    + r68;

% --- Reaction 69 ---
% O·⁻ + O₂ → O₃·⁻      k = 3.6×10^9 M⁻¹·s⁻¹
r69 = 3.6e9 * O_rad_minus * O2;
dO_rad_minus_dt = dO_rad_minus_dt - r69;
dO2_dt          = dO2_dt          - r69;
dO3_rad_dt      = dO3_rad_dt      + r69;

% --- Reaction 70 ---
% O·⁻ + H₂O → OH· + OH⁻      k = 1.8×10^5 M⁻¹·s⁻¹
% (If water is assumed constant [~55.5 M], include that factor)
r70 = 1.8e5 * O_rad_minus * 55.5;
dO_rad_minus_dt = dO_rad_minus_dt - r70;
dOH_rad_dt      = dOH_rad_dt      + r70;
dOH_minus_dt    = dOH_minus_dt    + r70;

% --- Reaction 71 ---
% O·⁻ + H₂O₂ → O₂·⁻ + H₂O      k = 4.0×10^8 M⁻¹·s⁻¹
r71 = 4.0e8 * O_rad_minus * H2O2;
dO_rad_minus_dt = dO_rad_minus_dt - r71;
dH2O2_dt        = dH2O2_dt        - r71;
dO2_super_dt    = dO2_super_dt    + r71;

% --- Reaction 72 ---
% (No reaction provided)

% --- Reaction 73 ---
% O₂·⁻ + H⁺ → HO₂·      k = 5×10^10 M⁻¹·s⁻¹
r73 = 5.0e10 * O2_super * H_plus;
dO2_super_dt = dO2_super_dt - r73;
dH_plus_dt   = dH_plus_dt   - r73;
dHO2_rad_dt  = dHO2_rad_dt  + r73;

% --- Reaction 74 ---
% (No reaction provided)

% --- Reaction 75 ---
% O₂·⁻ + OCl⁻ + H₂O → Cl· + 2OH⁻ + O₂      k = 2.0×10^8 M⁻¹·s⁻¹
r75 = 2.0e8 * O2_super * OCl_minus;
dO2_super_dt = dO2_super_dt - r75;
dOCl_minus_dt = dOCl_minus_dt - r75;
dCl_rad_dt   = dCl_rad_dt   + r75;
dOH_minus_dt = dOH_minus_dt + 2*r75;
dO2_dt       = dO2_dt       + r75;

% --- Reaction 76 ---
% O₃·⁻ → O·⁻ + O₂      k = 2.6×10^3 s⁻¹
r76 = 2.6e3 * O3_rad;
dO3_rad_dt     = dO3_rad_dt     - r76;
dO_rad_minus_dt = dO_rad_minus_dt + r76;
dO2_dt         = dO2_dt         + r76;

% --- Reaction 77 ---
% HO₂· + H₂O₂ → ·OH + H₂O + O₂      k = 3×10 M⁻¹·s⁻¹ | was mistakenly 3e10 M⁻¹·s⁻¹
r77 = 3.0 * HO2_rad * H2O2;
dHO2_rad_dt = dHO2_rad_dt - r77;
dH2O2_dt    = dH2O2_dt    - r77;
dOH_rad_dt  = dOH_rad_dt  + r77;
dO2_dt      = dO2_dt      + r77;

% --- Reaction 78 ---
% HO₂· + O₂·⁻ → HO₂⁻ + O₂      k = 9.7×10^7 M⁻¹·s⁻¹
r78 = 9.7e7 * HO2_rad * O2_super;
dHO2_rad_dt   = dHO2_rad_dt   - r78;
dO2_super_dt  = dO2_super_dt  - r78;
dHO2_minus_dt = dHO2_minus_dt + r78;
dO2_dt        = dO2_dt        + r78;

% --- Reaction 79 ---
% HO₂· → H⁺ + O₂·⁻      k = 7.0×10^5 s⁻¹
r79 = 7.0e5 * HO2_rad;
dHO2_rad_dt   = dHO2_rad_dt   - r79;
dH_plus_dt    = dH_plus_dt    + r79;
dO2_super_dt  = dO2_super_dt  + r79;

% --- Reaction 80 ---
% 2HO₂· + 2HOCl → 2Cl· + 2H₂O₂ + O₂      k = 7.5×10^6 M⁻¹·s⁻¹
r80 = 7.5e6 * HO2_rad * HOCl;
dHO2_rad_dt = dHO2_rad_dt - 2*r80;
dHOCl_dt    = dHOCl_dt    - 2*r80;
dCl_rad_dt  = dCl_rad_dt  + 2*r80;
dH2O2_dt    = dH2O2_dt    + 2*r80;
dO2_dt      = dO2_dt      + r80;

% --- Reaction 81 ---
% HO₂· + HO₂· → H₂O₂ + O₂      k = 8.3×10^5 M⁻¹·s⁻¹
r81 = 8.3e5 * (HO2_rad)^2;
dHO2_rad_dt = dHO2_rad_dt - 2*r81;
dH2O2_dt    = dH2O2_dt    + r81;
dO2_dt      = dO2_dt      + r81;

%% =======================================================
%% TABLE C3-2: Chlorine Photolysis Chemistry
%% =======================================================

% Reaction 47: OCl⁻ → O•⁻ + Cl•, ε = 66 M⁻¹cm⁻¹, φ = 0.97
r_C32_47 = (E0/L) * (1 - 10^(-a_total*L)) * (a_OCl / a_total) * phi_OCl_O_Cl;
dOCl_minus_dt = dOCl_minus_dt - r_C32_47;
dCl_rad_dt = dCl_rad_dt + r_C32_47;
dO_rad_minus_dt = dO_rad_minus_dt + r_C32_47;  % Now tracked

% Reaction 48: OCl⁻ → O(³P) + Cl⁻, φ = 0.074
r_C32_48 = (E0/L) * (1 - 10^(-a_total*L)) * (a_OCl / a_total) * phi_OCl_O3P;
dOCl_minus_dt = dOCl_minus_dt - r_C32_48;
dO3P_dt = dO3P_dt + r_C32_48;
dCl_minus_dt = dCl_minus_dt + r_C32_48;

% Reaction 49: OCl⁻ → O(¹D) + Cl⁻, φ = 0.133
r_C32_49 = (E0/L) * (1 - 10^(-a_total*L)) * (a_OCl / a_total) * phi_OCl_O1D;
dOCl_minus_dt = dOCl_minus_dt - r_C32_49;
dO1D_dt = dO1D_dt + r_C32_49;
dCl_minus_dt = dCl_minus_dt + r_C32_49;

%% =======================================================
%% TABLE C3-3:
%% =======================================================

% Reaction 54: ·OH + HO2· → O2 + H2O, k = 6.60e9 M⁻¹ s⁻¹
r_C33_54 = 6.60e9 * OH_rad * HO2_rad;
dOH_rad_dt = dOH_rad_dt - r_C33_54;
dHO2_rad_dt = dHO2_rad_dt - r_C33_54;
dO2_dt = dO2_dt + r_C33_54;

%% =======================================================
%% TABLE C3-6: ClOₓ Reactions
%% =======================================================

%% =======================================================
%% TABLE C3-6: ClOx Reactions (Including ClO·, Cl2O2, ClO2·, ClO2, ClO2-, ClO3-)
%% =======================================================

% -------------------------------
% Reaction 99: ClO· + NOM → Product 
r_C36_99 = 5.52e8 * ClO_rad * NOM;
dClO_rad_dt = dClO2_rad_dt - r_C36_99;


% Reaction 100: ClO· + ClO· → Cl2O2, k = 2.50e9 M⁻¹ s⁻¹
r_C36_100 = 2.50e9 * (ClO_rad)^2;
dClO_rad_dt = dClO_rad_dt - 2*r_C36_100;
dCl2O2_dt   = dCl2O2_dt   + r_C36_100;

% Reaction 101: ClO· + ClO2⁻ → ClO2· + OCl⁻, k = 9.40e8 M⁻¹ s⁻¹
r_C36_101 = 9.40e8 * ClO_rad * ClO2_minus;
dClO_rad_dt      = dClO_rad_dt      - r_C36_101;
dClO2_rad_dt     = dClO2_rad_dt     + r_C36_101;
dOCl_minus_dt    = dOCl_minus_dt    + r_C36_101;

% Reaction 102: ClO· + CO3²⁻ → OCl⁻ + CO3·⁻, k = 600 M⁻¹ s⁻¹
r_C36_102 = 600 * ClO_rad * CO3_2;
dClO_rad_dt      = dClO_rad_dt      - r_C36_102;
dOCl_minus_dt    = dOCl_minus_dt    + r_C36_102;
dCO3_rad_dt      = dCO3_rad_dt      + r_C36_102;

% Reaction 103: ClO· + ·OH → H⁺ + ClO2⁻, k = 1.00e9 M⁻¹ s⁻¹
r_C36_103 = 1.00e9 * ClO_rad * OH_rad;
dClO_rad_dt      = dClO_rad_dt      - r_C36_103;
dOH_rad_dt       = dOH_rad_dt       - r_C36_103;
dClO2_minus_dt   = dClO2_minus_dt   + r_C36_103;
dH_plus_dt       = dH_plus_dt       + r_C36_103;

% Reaction 104: 2ClO· + H2O → HOCl + H⁺ + ClO2⁻, k = 2.50e9 M⁻¹ s⁻¹
r_C36_104 = 2.50e9 * (ClO_rad)^2;
dClO_rad_dt      = dClO_rad_dt      - 2*r_C36_104;
dHOCl_dt         = dHOCl_dt         + r_C36_104;
dH_plus_dt       = dH_plus_dt       + r_C36_104;
dClO2_minus_dt   = dClO2_minus_dt   + r_C36_104;

% Reaction 105: 2ClO· + OH⁻ → OCl⁻ + H⁺ + ClO2⁻, k = 2.50e9 M⁻¹ s⁻¹
r_C36_105 = 2.50e9 * (ClO_rad)^2;
dClO_rad_dt      = dClO_rad_dt      - 2*r_C36_105;
dOCl_minus_dt    = dOCl_minus_dt    + r_C36_105;
dH_plus_dt       = dH_plus_dt       + r_C36_105;
dClO2_minus_dt   = dClO2_minus_dt   + r_C36_105;

% Reaction 106: Cl2O2 + H2O → ClO2⁻ + HOCl + H⁺, k = 1.00e4 s⁻¹
r_C36_106 = 1.00e4 * Cl2O2;
dCl2O2_dt      = dCl2O2_dt      - r_C36_106;
dClO2_minus_dt = dClO2_minus_dt + r_C36_106;
dHOCl_dt       = dHOCl_dt       + r_C36_106;
dH_plus_dt     = dH_plus_dt     + r_C36_106;

% Reaction 107: O(³P) + OCl⁻ → ClO2⁻, k = 9.40e9 M⁻¹ s⁻¹
r_C36_107 = 9.40e9 * O3P;  % (OCl⁻ is assumed implicitly in the rate)
dO3P_dt       = dO3P_dt       - r_C36_107;
dClO2_minus_dt = dClO2_minus_dt + r_C36_107;

% Reaction 108: ·OH + ClO2⁻ → ClO2· + OH⁻, k = 7.00e9 M⁻¹ s⁻¹
r_C36_108 = 7.00e9 * OH_rad * ClO2_minus;
dOH_rad_dt       = dOH_rad_dt       - r_C36_108;
dClO2_minus_dt   = dClO2_minus_dt   - r_C36_108;
dClO2_rad_dt     = dClO2_rad_dt     + r_C36_108;
dOH_minus_dt     = dOH_minus_dt     + r_C36_108;

% Reaction 109: Cl· + ClO2⁻ → ClO2· + Cl⁻, k = 7.00e9 M⁻¹ s⁻¹
r_C36_109 = 7.00e9 * Cl_rad * ClO2_minus;
dCl_rad_dt       = dCl_rad_dt       - r_C36_109;
dClO2_minus_dt   = dClO2_minus_dt   - r_C36_109;
dClO2_rad_dt     = dClO2_rad_dt     + r_C36_109;
dCl_minus_dt     = dCl_minus_dt     + r_C36_109;

% Reaction 110: Cl2·⁻ + ClO2⁻ → ClO2· + 2Cl⁻, k = 1.30e8 M⁻¹ s⁻¹
r_C36_110 = 1.30e8 * Cl2_rad_minus * ClO2_minus;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r_C36_110;
dClO2_minus_dt    = dClO2_minus_dt    - r_C36_110;
dClO2_rad_dt      = dClO2_rad_dt      + r_C36_110;
dCl_minus_dt      = dCl_minus_dt      + 2*r_C36_110;

% --- For Reaction 111 we now distinguish the atomic oxygen anion (O·⁻) from superoxide ---
% Reaction 111: O·⁻ + ClO2⁻ + H2O → 2OH⁻ + ClO2·, k = 1.95e8 M⁻¹ s⁻¹
r_C36_111 = 1.95e8 * O_rad_minus * ClO2_minus;
dO_rad_minus_dt = dO_rad_minus_dt - r_C36_111;
dClO2_minus_dt  = dClO2_minus_dt  - r_C36_111;
dOH_minus_dt    = dOH_minus_dt    + 2*r_C36_111;
dClO2_rad_dt   = dClO2_rad_dt   + r_C36_111;

% Reaction 112: O2·⁻ + ClO2⁻ → Product, k = 40 M⁻¹ s⁻¹
r_C36_112 = 40 * O2_super * ClO2_minus;
dO2_super_dt   = dO2_super_dt   - r_C36_112;
dClO2_minus_dt = dClO2_minus_dt - r_C36_112;

% Reaction 113: ClO2⁻ + HOCl + H⁺ → Cl2O2 + H2O, k = 2.86e6 M⁻¹ s⁻¹
r_C36_113 = 2.86e6 * ClO2_minus * HOCl * H_plus;
dClO2_minus_dt = dClO2_minus_dt - r_C36_113;
dHOCl_dt       = dHOCl_dt       - r_C36_113;
dH_plus_dt     = dH_plus_dt     - r_C36_113;
dCl2O2_dt     = dCl2O2_dt      + r_C36_113;

% Reaction 114: Cl2O2 + ClO2⁻ → Cl⁻ + 2ClO2, k = 8.10e5 M⁻¹ s⁻¹
r_C36_114 = 8.10e5 * Cl2O2 * ClO2_minus;
dCl2O2_dt      = dCl2O2_dt      - r_C36_114;
dClO2_minus_dt = dClO2_minus_dt - r_C36_114;
dCl_minus_dt   = dCl_minus_dt   + r_C36_114;
dClO2_dt       = dClO2_dt       + 2*r_C36_114;

% Reaction 115: 2HOCl + ClO2⁻ → ClO3⁻ + Cl2 + H2O, k = 2.10e3 M⁻¹ s⁻¹
r_C36_115 = 2.10e3 * HOCl * ClO2_minus;
dHOCl_dt       = dHOCl_dt       - 2*r_C36_115; % two HOCl are consumed
dClO2_minus_dt = dClO2_minus_dt - r_C36_115;
dClO3_minus_dt = dClO3_minus_dt + r_C36_115;
dCl2_dt        = dCl2_dt        + r_C36_115;

% Reaction 116: Cl2 + ClO2⁻ → Cl2O2 + Cl⁻, k = 1.61e6 M⁻¹ s⁻¹
r_C36_116 = 1.61e6 * Cl2 * ClO2_minus;
dCl2_dt        = dCl2_dt        - r_C36_116;
dClO2_minus_dt = dClO2_minus_dt - r_C36_116;
dCl2O2_dt      = dCl2O2_dt      + r_C36_116;
dCl_minus_dt   = dCl_minus_dt   + r_C36_116;

% Reaction 117: CO3·⁻ + ClO2⁻ → ClO2· + CO3²⁻, k = 3.40e7 M⁻¹ s⁻¹
r_C36_117 = 3.40e7 * CO3_rad * ClO2_minus;
dCO3_rad_dt    = dCO3_rad_dt    - r_C36_117;
dClO2_minus_dt = dClO2_minus_dt - r_C36_117;
dClO2_rad_dt   = dClO2_rad_dt   + r_C36_117;
dCO3_2_dt      = dCO3_2_dt      + r_C36_117;

% Reaction 118: ·OH + ClO2· → ClO3⁻ + H⁺, k = 4.00e9 M⁻¹ s⁻¹
r_C36_118 = 4.00e9 * OH_rad * ClO2_rad;
dOH_rad_dt       = dOH_rad_dt       - r_C36_118;
dClO2_rad_dt     = dClO2_rad_dt     - r_C36_118;
dClO3_minus_dt   = dClO3_minus_dt   + r_C36_118;
dH_plus_dt       = dH_plus_dt       + r_C36_118;

% Reaction 119: O·⁻ + ClO2· → ClO3⁻, k = 2.70e9 M⁻¹ s⁻¹
r_C36_119 = 2.70e9 * O_rad_minus * ClO2_rad;
dO_rad_minus_dt = dO_rad_minus_dt - r_C36_119;
dClO2_rad_dt    = dClO2_rad_dt    - r_C36_119;
dClO3_minus_dt  = dClO3_minus_dt  + r_C36_119;

% Reaction 120: Cl· + ClO2· → Product, k = 4.00e9 M⁻¹ s⁻¹
r_C36_120 = 4.00e9 * Cl_rad * ClO2_rad;
dCl_rad_dt     = dCl_rad_dt     - r_C36_120;
dClO2_rad_dt   = dClO2_rad_dt   - r_C36_120;
% (Product not tracked)

% Reaction 121: ClO2· + ClO· + H2O → ClO3⁻ + HOCl + H⁺, k = 9.40e9 M⁻¹ s⁻¹
r_C36_121 = 9.40e9 * ClO2_rad * ClO_rad;
dClO2_rad_dt   = dClO2_rad_dt   - r_C36_121;
dClO_rad_dt    = dClO_rad_dt    - r_C36_121;
dClO3_minus_dt = dClO3_minus_dt + r_C36_121;
dHOCl_dt       = dHOCl_dt       + r_C36_121;
dH_plus_dt     = dH_plus_dt     + r_C36_121;

% Reaction 122: O2·⁻ + ClO2· → ClO2⁻ + O2, k = 3.15e9 M⁻¹ s⁻¹
r_C36_122 = 3.15e9 * O2_super * ClO2_rad;
dO2_super_dt   = dO2_super_dt   - r_C36_122;
dClO2_rad_dt   = dClO2_rad_dt   - r_C36_122;
dClO2_minus_dt = dClO2_minus_dt + r_C36_122;
% (O2 is not explicitly tracked)

% Reaction 123: HO2⁻ + ClO2· → ClO2⁻ + HO2·, k = 9.57e4 M⁻¹ s⁻¹
r_C36_123 = 9.57e4 * HO2_minus * ClO2_rad;
dHO2_minus_dt  = dHO2_minus_dt - r_C36_123;
dClO2_rad_dt   = dClO2_rad_dt   - r_C36_123;
dClO2_minus_dt = dClO2_minus_dt + r_C36_123;
dHO2_rad_dt    = dHO2_rad_dt    + r_C36_123;

% Reaction 124: HO2· + ClO2⁻ → Product, k = 1.00e6 M⁻¹ s⁻¹
r_C36_124 = 1.00e6 * HO2_rad * ClO2_minus;
dHO2_rad_dt    = dHO2_rad_dt    - r_C36_124;
dClO2_minus_dt = dClO2_minus_dt - r_C36_124;
% (Product not tracked)

% Reaction 125: Cl· + ClO2⁻ → Cl2O2, k = 7.80e9 M⁻¹ s⁻¹
r_C36_125 = 7.80e9 * Cl_rad * ClO2_minus;
dCl_rad_dt     = dCl_rad_dt     - r_C36_125;
dClO2_minus_dt = dClO2_minus_dt - r_C36_125;
dCl2O2_dt      = dCl2O2_dt      + r_C36_125;

% Reaction 126: Cl2·⁻ + ClO2· → Product, k = 1.00e9 M⁻¹ s⁻¹
r_C36_126 = 1.00e9 * Cl2_rad_minus * ClO2_rad;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r_C36_126;
dClO2_rad_dt      = dClO2_rad_dt      - r_C36_126;
% (Product not tracked)

% Reaction 127: O3·⁻ + ClO2⁻ → ClO3⁻ + O2, k = 1.80e5 M⁻¹ s⁻¹
r_C36_127 = 1.80e5 * O3_rad_minus * ClO2_minus;
dO3_rad_minus_dt = dO3_rad_minus_dt - r_C36_127;
dClO2_minus_dt   = dClO2_minus_dt   - r_C36_127;
dClO3_minus_dt   = dClO3_minus_dt   + r_C36_127;
% (O2 is not explicitly tracked)

% Reaction 128: O(³P) + ClO2⁻ → ClO3⁻, k = 1.59e10 M⁻¹ s⁻¹
r_C36_128 = 1.59e10 * O3P * ClO2_minus;
dO3P_dt        = dO3P_dt        - r_C36_128;
dClO2_minus_dt = dClO2_minus_dt - r_C36_128;
dClO3_minus_dt = dClO3_minus_dt + r_C36_128;

% Reaction 129: ·OH + ClO3⁻ → Product, k = 1.00e6 M⁻¹ s⁻¹
r_C36_129 = 1.00e6 * OH_rad * ClO3_minus;
dOH_rad_dt     = dOH_rad_dt     - r_C36_129;
dClO3_minus_dt = dClO3_minus_dt - r_C36_129;
% (Product not tracked)

% Reaction 130: Cl· + ClO3⁻ → Product, k = 1.00e6 M⁻¹ s⁻¹
r_C36_130 = 1.00e6 * Cl_rad * ClO3_minus;
dCl_rad_dt     = dCl_rad_dt     - r_C36_130;
dClO3_minus_dt = dClO3_minus_dt - r_C36_130;
% (Product not tracked)

% Reaction 131: O2·⁻ + ClO3⁻ → Product, k = 0.003 M⁻¹ s⁻¹
r_C36_131 = 0.003 * O2_super * ClO3_minus;
dO2_super_dt   = dO2_super_dt   - r_C36_131;
dClO3_minus_dt = dClO3_minus_dt - r_C36_131;
% (Product not tracked)


%% =======================================================
%% TABLE C3-7: ClOH•⁻ Reactions
%% =======================================================

% Reaction 132: ClOH•⁻ + H⁺ → Cl• + H2O, k = 2.10e10 M⁻¹ s⁻¹ (duplicate)

% Reaction 133:  ClOH•⁻ + H⁺ → Cl• + H2O, k = 2.10×10^10 M⁻¹ s⁻¹
r_C37_133 = 2.10e10 * ClOH_rad_minus * H_plus;
dClOH_rad_minus_dt = dClOH_rad_minus_dt - r_C37_133;
dH_plus_dt = dH_plus_dt - r_C37_133;
dCl_rad_dt = dCl_rad_dt + r_C37_133;

% Reaction 135:  ClOH•⁻ → Cl• + OH⁻, k = 15.1 s⁻¹
r_C37_135 = 15.1 * ClOH_rad_minus;
dClOH_rad_minus_dt = dClOH_rad_minus_dt - r_C37_135;
dCl_rad_dt = dCl_rad_dt + r_C37_135;
dOH_rad_dt = dOH_rad_dt + r_C37_135;

%% =======================================================
%% TABLE C3-8: HO2· and O2·⁻ Reactions
%% =======================================================

% Reaction 137: HO2· + H2O2 → O2 + ·OH + H2O, k = 3.00 M⁻¹ s⁻¹
r_C38_137 = 3.00 * H2O2 * HO2_rad;
dH2O2_dt = dH2O2_dt - r_C38_137;
dHO2_rad_dt = dHO2_rad_dt - r_C38_137;
dOH_rad_dt = dOH_rad_dt + r_C38_137;
dO2_dt = dO2_dt + r_C38_137;
% Reaction 138: HO2· + O2·⁻ → O2 + HO2⁻, k = 9.70e7 M⁻¹ s⁻¹
r_C38_138 = 9.70e7 * HO2_rad * O2_super;
dHO2_rad_dt = dHO2_rad_dt - r_C38_138;
dO2_super_dt = dO2_super_dt - r_C38_138;
dHO2_minus_dt = dHO2_minus_dt + r_C38_138;
dO2_dt = dO2_dt + r_C38_138;
% Reaction 139: HO2· + HOCl → Cl• + H2O + O2, k = 7.50e6 M⁻¹ s⁻¹
r_C38_139 = 7.50e6 * HO2_rad * HOCl;
dHO2_rad_dt = dHO2_rad_dt - r_C38_139;
dHOCl_dt = dHOCl_dt - r_C38_139;
dCl_rad_dt = dCl_rad_dt + r_C38_139;
dO2_dt = dO2_dt + r_C38_139;
% Reaction 140: O2·⁻ + OCl⁻ + H2O → Cl• + 2OH⁻ + O2, k = 2.00e8 M⁻¹ s⁻¹
r_C38_140 = 2.00e8 * O2_super * OCl_minus;
dO2_super_dt = dO2_super_dt - r_C38_140;
dOCl_minus_dt = dOCl_minus_dt - r_C38_140;
dCl_rad_dt = dCl_rad_dt + r_C38_140;
dOH_minus_dt = dOH_minus_dt + 2*r_C38_140;
dO2_dt = dO2_dt + r_C38_140;
% Reaction 141: O2·⁻ + H2O2 → O2 + ·OH + OH⁻, k = 1.30e-1 M⁻¹ s⁻¹
r_C38_141 = 1.30e-1 * O2_super * H2O2;
dO2_super_dt = dO2_super_dt - r_C38_141;
dH2O2_dt = dH2O2_dt - r_C38_141;
dOH_rad_dt = dOH_rad_dt + r_C38_141;
dO2_dt = dO2_dt + r_C38_141;
dOH_minus_dt = dOH_minus_dt + r_C38_141;
% Reaction 142: O2·⁻ + HOCl → Cl• + OH⁻ + O2, k = 7.50e6 M⁻¹ s⁻¹
r_C38_142 = 7.50e6 * O2_super * HOCl;
dO2_super_dt = dO2_super_dt - r_C38_142;
dHOCl_dt = dHOCl_dt - r_C38_142;
dCl_rad_dt = dCl_rad_dt + r_C38_142;
dOH_minus_dt = dOH_minus_dt + r_C38_142;
dO2_dt = dO2_dt + r_C38_142;
% Reaction 143: O2·⁻ + CO3·⁻ → CO3²⁻ + O2, k = 6.00e8 M⁻¹ s⁻¹
r_C38_143 = 6.00e8 * O2_super * CO3_rad;
dO2_super_dt = dO2_super_dt - r_C38_143;
dCO3_rad_dt = dCO3_rad_dt - r_C38_143;
dCO3_2_dt = dCO3_2_dt + r_C38_143;
dO2_dt = dO2_dt + r_C38_143;
% Reaction 144: O2·⁻ + Cl⁻ → Product, k = 140 M⁻¹ s⁻¹
r_C38_144 = 140 * O2_super * Cl_minus;
dO2_super_dt = dO2_super_dt - r_C38_144;
% Reaction 145: (Not defined; assume 0)
r_C38_145 = 0;

%% =======================================================
%% TABLE C3-9: CO3·⁻ Reactions
%% =======================================================

% Reaction 146: CO3·⁻ + H2O2 → HCO3⁻ + HO2·, k = 4.30e5 M⁻¹ s⁻¹
r_C39_146 = 4.30e5 * CO3_rad * H2O2;
dCO3_rad_dt = dCO3_rad_dt - r_C39_146;
dHCO3_minus_dt = dHCO3_minus_dt + r_C39_146;
dHO2_rad_dt = dHO2_rad_dt + r_C39_146;
% Reaction 147: CO3·⁻ + HO2⁻ → CO3²⁻ + HO2·, k = 3.00e7 M⁻¹ s⁻¹
r_C39_147 = 3.00e7 * CO3_rad * HO2_minus;
dCO3_rad_dt = dCO3_rad_dt - r_C39_147;
dCO3_2_dt = dCO3_2_dt + r_C39_147;
dHO2_rad_dt = dHO2_rad_dt + r_C39_147;
% Reaction 148: CO3·⁻ + OCl⁻ → ClO· + CO3²⁻, k = 5.70e5 M⁻¹ s⁻¹
r_C39_148 = 5.70e5 * CO3_rad * OCl_minus;
dCO3_rad_dt = dCO3_rad_dt - r_C39_148;
dClO_rad_dt = dClO_rad_dt + r_C39_148;
dCO3_2_dt = dCO3_2_dt + r_C39_148;
% Reaction 149: CO3·⁻ + NOM → Product, k = 6.96e5 M⁻¹ s⁻¹ (NOM not tracked)
r_C39_149 = 6.96e5 * CO3_rad * NOM;
dCO3_rad_dt = dCO3_rad_dt - r_C39_149;
% Reaction 150: (Not defined; assume 0)
r_C39_150 = 0;

%% =======================================================
%% TABLE C3-10: Equilibrium Reactions (Phosphate/Carbonate)
%% =======================================================

% Reaction 151: H2O2 → H⁺ + HO2⁻, k_f = 1.26e-1 s⁻¹
r_C310_151 = 1.26e-1 * H2O2;
dH2O2_dt = dH2O2_dt - r_C310_151;
dH_plus_dt = dH_plus_dt + r_C310_151;
dHO2_minus_dt = dHO2_minus_dt + r_C310_151;
% Reaction 152: H⁺ + HO2⁻ → H2O2, k = 5.00e10 M⁻¹ s⁻¹
r_C310_152 = 5.00e10 * H_plus * HO2_minus;
dH_plus_dt = dH_plus_dt - r_C310_152;
dHO2_minus_dt = dHO2_minus_dt - r_C310_152;
dH2O2_dt = dH2O2_dt + r_C310_152;
% Reactions 153-157 (Phosphate reactions): not tracked; set to 0.
r_C310_153 = 0; r_C310_154 = 0; r_C310_155 = 0; r_C310_156 = 0; r_C310_157 = 0;
% Reaction 164: HCO3⁻ + H⁺ → H2CO3, k = 5.00e10 M⁻¹ s⁻¹
r_C310_164 = 5.00e10 * HCO3_minus * H_plus;
dHCO3_minus_dt = dHCO3_minus_dt - r_C310_164;
dH_plus_dt = dH_plus_dt - r_C310_164;
dH2CO3_dt = dH2CO3_dt + r_C310_164;
% Reaction 165: H2CO3 → HCO3⁻ + H⁺, k = 5.00e5 s⁻¹
r_C310_165 = 5.00e5 * H2CO3;
dH2CO3_dt = dH2CO3_dt - r_C310_165;
dHCO3_minus_dt = dHCO3_minus_dt + r_C310_165;
dH_plus_dt = dH_plus_dt + r_C310_165;
% Reaction 166: CO3²⁻ + H⁺ → HCO3⁻, k = 5.00e10 M⁻¹ s⁻¹
r_C310_166 = 5.00e10 * CO3_2 * H_plus;
dCO3_2_dt = dCO3_2_dt - r_C310_166;
dH_plus_dt = dH_plus_dt - r_C310_166;
dHCO3_minus_dt = dHCO3_minus_dt + r_C310_166;
% Reaction 167: HCO3⁻ → CO3²⁻ + H⁺, k = 2.50 s⁻¹
r_C310_167 = 2.50 * HCO3_minus;
dHCO3_minus_dt = dHCO3_minus_dt - r_C310_167;
dCO3_2_dt = dCO3_2_dt + r_C310_167;
dH_plus_dt = dH_plus_dt + r_C310_167;


%% =======================================================
%% TABLE C3-11: O₃ and O₃·⁻ Reactions
%% =======================================================

% Reaction 169: O(³P) + O₂ → O₃, k = 4.00e9 M⁻¹ s⁻¹
r_C311_169 = 4.00e9 * O3P * O2;
dO3P_dt = dO3P_dt - r_C311_169;
dO2_dt = dO2_dt - r_C311_169;
dO3_var_dt = dO3_var_dt + r_C311_169;
% Reaction 170: O₃ → O(³P) + O₂, k = 4.50e-6 s⁻¹
r_C311_170 = 4.50e-6 * O3_var;
dO3_var_dt = dO3_var_dt - r_C311_170;
dO3P_dt = dO3P_dt + r_C311_170;
dO2_dt = dO2_dt + r_C311_170;
% Reaction 171: O₃ + ClO2⁻ → O₃·⁻ + ClO2·, k = 4.00e6 M⁻¹ s⁻¹
r_C311_171 = 4.00e6 * O3_var * ClO2_minus;
dO3_var_dt = dO3_var_dt - r_C311_171;
dClO2_minus_dt = dClO2_minus_dt - r_C311_171;
dO3_rad_dt = dO3_rad_dt + r_C311_171;
dClO2_rad_dt = dClO2_rad_dt + r_C311_171;
% Reaction 172: O₃ + ClO2· → O₂ + ClO3⁻, k = 1.23e3 M⁻¹ s⁻¹
r_C311_172 = 1.23e3 * O3_var * ClO2_rad;
dO3_var_dt = dO3_var_dt - r_C311_172;
dClO2_rad_dt = dClO2_rad_dt - r_C311_172;
dO2_dt = dO2_dt + r_C311_172;
dClO3_minus_dt = dClO3_minus_dt + r_C311_172;
% Reaction 173: O₃ + OCl⁻ → 2O₂ + Cl⁻, k = 1.10e2 M⁻¹ s⁻¹
r_C311_173 = 1.10e2 * O3_var * OCl_minus;
dO3_var_dt = dO3_var_dt - r_C311_173;
dOCl_minus_dt = dOCl_minus_dt - r_C311_173;
dO2_dt = dO2_dt + 2*r_C311_173;
dCl_minus_dt = dCl_minus_dt + r_C311_173;
% Reaction 174: O₃ + OCl⁻ → O₂ + ClO2⁻, k = 30 M⁻¹ s⁻¹
r_C311_174 = 30 * O3_var * OCl_minus;
dO3_var_dt = dO3_var_dt - r_C311_174;
dOCl_minus_dt = dOCl_minus_dt - r_C311_174;
dO2_dt = dO2_dt + r_C311_174;
dClO2_minus_dt = dClO2_minus_dt + r_C311_174;
% Reaction 175: O₃ + Cl⁻ → O₂ + OCl⁻, k = 1.60e-3 M⁻¹ s⁻¹
r_C311_175 = 1.60e-3 * O3_var * Cl_minus;
dO3_var_dt = dO3_var_dt - r_C311_175;
dCl_minus_dt = dCl_minus_dt - r_C311_175;
dO2_dt = dO2_dt + r_C311_175;
dOCl_minus_dt = dOCl_minus_dt + r_C311_175;
% Reaction 176: O₃ + Cl2•⁻ → Product, k = 9.00e7 M⁻¹ s⁻¹
r_C311_176 = 9.00e7 * O3_var * Cl2_rad_minus;
dO3_var_dt = dO3_var_dt - r_C311_176;
dCl2_rad_minus_dt = dCl2_rad_minus_dt - r_C311_176;
% Reaction 177: O₃ + ·OH → O₂ + HO2·, k = 1.10e8 M⁻¹ s⁻¹
r_C311_177 = 1.10e8 * O3_var * OH_rad;
dO3_var_dt = dO3_var_dt - r_C311_177;
dOH_rad_dt = dOH_rad_dt - r_C311_177;
dO2_dt = dO2_dt + r_C311_177;
dHO2_rad_dt = dHO2_rad_dt + r_C311_177;
% Reaction 178: O₃ + OH⁻ → O₂ + HO2⁻, k = 14.2 M⁻¹ s⁻¹
r_C311_178 = 14.2 * O3_var * OH_minus;
dO3_var_dt = dO3_var_dt - r_C311_178;
dOH_minus_dt = dOH_minus_dt - r_C311_178;
dO2_dt = dO2_dt + r_C311_178;
dHO2_minus_dt = dHO2_minus_dt + r_C311_178;
% Reaction 179: O₃ + H· → O₂ + ·OH, k = 2.20e10 M⁻¹ s⁻¹
r_C311_179 = 2.20e10 * O3_var * H_atom;
dO3_var_dt = dO3_var_dt - r_C311_179;
dH_atom_dt = dH_atom_dt - r_C311_179;
dO2_dt = dO2_dt + r_C311_179;
dOH_rad_dt = dOH_rad_dt + r_C311_179;
% Reaction 180: O₃ + HO2· → O₂ + H⁺ + O₃·⁻, k = 1.60e9 M⁻¹ s⁻¹
r_C311_180 = 1.60e9 * O3_var * HO2_rad;
dO3_var_dt = dO3_var_dt - r_C311_180;
dHO2_rad_dt = dHO2_rad_dt - r_C311_180;
dO3_rad_dt = dO3_rad_dt + r_C311_180;
dH_plus_dt = dH_plus_dt + r_C311_180;
dO2_dt = dO2_dt + r_C311_180;
% Reaction 181: O₃ + HO2⁻ → O₃·⁻ + HO2·, k = 5.50e6 M⁻¹ s⁻¹
r_C311_181 = 5.50e6 * O3_var * HO2_minus;
dO3_var_dt = dO3_var_dt - r_C311_181;
dHO2_minus_dt = dHO2_minus_dt - r_C311_181;
dO3_rad_dt = dO3_rad_dt + r_C311_181;
dHO2_rad_dt = dHO2_rad_dt + r_C311_181;
% Reaction 182: O₃·⁻ + ·OH → O₂·⁻ + HO2·, k = 8.50e9 M⁻¹ s⁻¹
r_C311_182 = 8.50e9 * O3_rad * OH_rad;
dO3_rad_dt = dO3_rad_dt - r_C311_182;
dOH_rad_dt = dOH_rad_dt - r_C311_182;
dO2_super_dt = dO2_super_dt + r_C311_182;
dHO2_rad_dt = dHO2_rad_dt + r_C311_182;
% Reaction 183: O₃·⁻ + H⁺ → O₂ + ·OH, k = 5.20e10 M⁻¹ s⁻¹
r_C311_183 = 5.20e10 * O3_rad * H_plus;
dO3_rad_dt = dO3_rad_dt - r_C311_183;
dH_plus_dt = dH_plus_dt - r_C311_183;
dO2_dt = dO2_dt + r_C311_183;
dOH_rad_dt = dOH_rad_dt + r_C311_183;
% Reaction 184: O₃·⁻ + ClO· → O₃ + OCl⁻, k = 1.00e9 M⁻¹ s⁻¹
r_C311_184 = 1.00e9 * O3_rad * ClO_rad;
dO3_rad_dt = dO3_rad_dt - r_C311_184;
dClO_rad_dt = dClO_rad_dt - r_C311_184;
dO3_var_dt = dO3_var_dt + r_C311_184;
dOCl_minus_dt = dOCl_minus_dt + r_C311_184;
% Reaction 185: O₃·⁻ → O₂ + O·⁻, k = 4.28e3 s⁻¹
r_C311_185 = 4.28e3 * O3_rad;
dO3_rad_dt = dO3_rad_dt - r_C311_185;
dO2_dt = dO2_dt + r_C311_185;
dO2_super_dt = dO2_super_dt + r_C311_185;
% Reaction 186: O₃·⁻ + O₃·⁻ → Product, k = 9.00e8 M⁻¹ s⁻¹
r_C311_186 = 9.00e8 * (O3_rad)^2;
dO3_rad_dt = dO3_rad_dt - r_C311_186;
% Reaction 187: O₃·⁻ + O·⁻ → 2O₂·⁻, k = 7.00e8 M⁻¹ s⁻¹
r_C311_187 = 7.00e8 * O3_rad * O2_super;
dO3_rad_dt = dO3_rad_dt - r_C311_187;
dO2_super_dt = dO2_super_dt + r_C311_187;
% Reaction 188: O₂·⁻ + O₃ → O₃·⁻ + O₂, k = 1.60e9 M⁻¹ s⁻¹
r_C311_188 = 1.60e9 * O2_super * O3_var;
dO2_super_dt = dO2_super_dt - r_C311_188;
dO3_rad_dt = dO3_rad_dt + r_C311_188;
dO2_dt = dO2_dt + r_C311_188;

%% =======================================================
%% TABLE C3-12: Cl₂ and Cl₃⁻ Reactions
%% =======================================================

% Reaction 189: Cl₂ + H2O2 → 2HCl + O₂, k = 1.30e4 M⁻¹ s⁻¹
r_C312_189 = 1.30e4 * Cl2 * H2O2;
dCl2_dt = dCl2_dt - r_C312_189;
dHCl_dt = dHCl_dt + 2*r_C312_189;
dO2_dt = dO2_dt + r_C312_189;

% Reaction 190: Cl₂ + O₂·⁻ → Cl2•⁻ + O₂, k = 1.00e9 M⁻¹ s⁻¹
r_C312_190 = 1.00e9 * Cl2 * O2_super;
dCl2_dt = dCl2_dt - r_C312_190;
dO2_super_dt = dO2_super_dt - r_C312_190;
dCl2_rad_minus_dt = dCl2_rad_minus_dt + r_C312_190;
dO2_dt = dO2_dt + r_C312_190;

% Reaction 191: Cl₂ + HO2· → Cl2•⁻ + H⁺ + O₂, k = 1.00e9 M⁻¹ s⁻¹
r_C312_191 = 1.00e9 * Cl2 * HO2_rad;
dCl2_dt = dCl2_dt - r_C312_191;
dHO2_rad_dt = dHO2_rad_dt - r_C312_191;
dH_plus_dt = dH_plus_dt + r_C312_191;
dCl2_rad_minus_dt = dCl2_rad_minus_dt + r_C312_191;
% dO2_dt = dO2_dt + r_C312_191;

% Reaction 192: Cl₃⁻ + HO2· → Cl2•⁻ + H⁺ + Cl⁻ + O₂, k = 1.00e9 M⁻¹ s⁻¹
r_C312_192 = 1.00e9 * Cl3_minus * HO2_rad;
dCl3_minus_dt = dCl3_minus_dt - r_C312_192;
dHO2_rad_dt = dHO2_rad_dt - r_C312_192;
dH_plus_dt = dH_plus_dt + r_C312_192;
dCl_minus_dt = dCl_minus_dt + r_C312_192;
dCl2_rad_minus_dt = dCl2_rad_minus_dt + r_C312_192;
dO2_dt = dO2_dt + r_C312_192;

% Reaction 193: Cl₃⁻ + O₂·⁻ → Cl2•⁻ + Cl⁻ + O₂, k = 3.80e9 M⁻¹ s⁻¹
r_C312_193 = 3.80e9 * Cl3_minus * O2_super;
dCl3_minus_dt = dCl3_minus_dt - r_C312_193;
dO2_super_dt = dO2_super_dt - r_C312_193;
dCl_minus_dt = dCl_minus_dt + r_C312_193;
dCl2_rad_minus_dt = dCl2_rad_minus_dt + r_C312_193;
dO2_dt = dO2_dt + r_C312_193;

% Reaction 194: Cl₂ + Cl• → Cl₃• (assumed to rapidly dissociate; no net change)
r_C312_194 = 5.30e8 * Cl2 * Cl_rad;
dCl2_dt = dCl2_dt - r_C312_194;

% Reaction 195: Cl₂ + OH⁻ → HOCl + Cl⁻, k = 1.00e9 M⁻¹ s⁻¹
r_C312_195 = 1.00e9 * Cl2 * OH_minus;
dCl2_dt = dCl2_dt - r_C312_195;
dOH_minus_dt = dOH_minus_dt - r_C312_195;
dHOCl_dt = dHOCl_dt + r_C312_195;
dCl_minus_dt = dCl_minus_dt + r_C312_195;

% Reactions 196–197: (Not provided; assume 0)
r_C312_196 = 0;
r_C312_197 = 0;

%% =======================================================
%% TABLE C3-13: Miscellaneous Reactions
%% =======================================================

% Reaction 198: HClOH• → ClOH•⁻ + H⁺, k = 1.00e8 s⁻¹
r_C313_198 = 1.00e8 * HClOH_rad;
dHClOH_rad_dt = dHClOH_rad_dt - r_C313_198;
dClOH_rad_minus_dt = dClOH_rad_minus_dt + r_C313_198;
dH_plus_dt = dH_plus_dt + r_C313_198;

% Reaction 199: HClOH• → Cl• + H2O, k = 1.00e2 s⁻¹
r_C313_199 = 1.00e2 * HClOH_rad;
dHClOH_rad_dt = dHClOH_rad_dt - r_C313_199;
dCl_rad_dt = dCl_rad_dt + r_C313_199;

% Reaction 200: HClOH• + Cl⁻ → Cl2•⁻ + H2O, k = 5.00e9 M⁻¹ s⁻¹
r_C313_200 = 5.00e9 * HClOH_rad * Cl_minus;
dHClOH_rad_dt = dHClOH_rad_dt - r_C313_200;
dCl_minus_dt = dCl_minus_dt - r_C313_200;
dCl2_rad_minus_dt = dCl2_rad_minus_dt + r_C313_200;

% Reaction 201: HOCl + H2O2 → H⁺ + Cl⁻ + H2O + O2, k = 1.10e4 M⁻¹ s⁻¹
r_C313_201 = 1.10e4 * HOCl * H2O2;
dHOCl_dt = dHOCl_dt - r_C313_201;
dH2O2_dt = dH2O2_dt - r_C313_201;
dH_plus_dt = dH_plus_dt + r_C313_201;
dCl_minus_dt = dCl_minus_dt + r_C313_201;
dO2_dt = dO2_dt + r_C313_201;

% Reaction 202: OCl⁻ + H2O2 → Cl⁻ + H2O + O2, k = 1.70e5 M⁻¹ s⁻¹
r_C313_202 = 1.70e5 * OCl_minus * H2O2;
dOCl_minus_dt = dOCl_minus_dt - r_C313_202;
dH2O2_dt = dH2O2_dt - r_C313_202;
dCl_minus_dt = dCl_minus_dt + r_C313_202;
dO2_dt = dO2_dt + r_C313_202;

% Reaction 203: (Not provided; assume 0)
r_C313_203 = 0;

% % Reaction 204: Cl⁻ + HOCl + H⁺ → Cl₂ + H2O, k = 0.182 M⁻¹ s⁻¹
% r_C313_204 = 0.182 * H_plus * HOCl * Cl_minus;
% dH_plus_dt = dH_plus_dt - r_C313_204;
% dHOCl_dt = dHOCl_dt - r_C313_204;
% dCl_minus_dt = dCl_minus_dt - r_C313_204;
% dCl2_dt = dCl2_dt + r_C313_204;

% Reaction 205: (Not provided; assume 0)
r_C313_205 = 0;

% Reaction from Chen's model: OH radical background scagenging
r_background_scavenging = app.sc * OH_rad;
dOH_rad_dt = dOH_rad_dt - r_background_scavenging;


% OH exposure
dOHexposure_dt = OH_rad;  % Integrating OH_rad over time

%% (End of report's model tables)
%% =======================================================




%% Assemble the final derivative vector (52×1)
dydt = [ dNH2Cl_dt;        % 1
         dNHCl2_dt;        % 2
         dHOCl_dt;         % 3
         dOCl_minus_dt;    % 4
         dNHCl_rad_dt;     % 5
         dNH2_rad_dt;      % 6
         dNCl2_rad_dt;     % 7
         dNCl3_dt;         % 8
         dNH3_dt;          % 9
         dNH4_dt;          % 10
         dH2O2_dt;         % 11
         dOH_rad_dt;       % 12
         dHO2_rad_dt;      % 13
         dHO2_minus_dt;    % 14
         dO3P_dt;          % 15
         dO1D_dt;          % 16
         dCl_rad_dt;       % 17
         dCl2_rad_minus_dt;% 18
         dCl2_dt;          % 19
         dClO_rad_dt;      % 20
         dClO2_minus_dt;   % 21
         dClO2_rad_dt;     % 22
         dClO3_minus_dt;   % 23
         dCl2O2_dt;        % 24
         dClOH_rad_minus_dt;% 25
         dCl3_minus_dt;    % 26
         dNOH_dt;          % 27
         dNO_rad_dt;       % 28
         dNO2_rad_dt;      % 29
         dNO2_minus_dt;    % 30
         dN2O3_dt;         % 31
         dN2O4_dt;         % 32
         dNO3_minus_dt;    % 33
         ddioxane_dt;      % 34
         dDOC_dt;          % 35
         dO2_dt;           % 36
         dH_plus_dt;       % 37
         dOH_minus_dt;     % 38
         dHCO3_minus_dt;   % 39
         dH2CO3_dt;        % 40
         dCO3_2_dt;        % 41
         dCl_minus_dt;     % 42
         dHCl_dt;          % 43
         dNH2O2_rad_dt;    % 44
         dNHClO2_rad_dt;   % 45
         dN2O_dt;          % 46
         dCO3_rad_dt;      % 47
         dN2_dt;           % 48
         dO_rad_minus_dt;  % 49
         dHPO4_2_dt;       % 50
         dH2PO4_minus_dt;  % 51
         dHPO4_rad_dt;
         dHClOH_rad_dt;
         dO2_super_dt;
         dO3_var_dt;
         dO3_rad_dt;
         dH_atom_dt; %H atom
         dClO2_dt; %ClO2
         dO3_rad_minus_dt; %O3 radical
         dNOM_dt;
         dcaffeine_dt; %caffeine
         dsucralose_dt;
         dOHexposure_dt; %OH exposure
         ];   % 63
         % Clamp the pH
dydt(37) = 0;  % Fix [H⁺] to its initial value so pH remains constant
dydt(38) = 0;  % Fix [OH-] to its initial value so pH remains constant

        %% Clamp alkalinity (Fix HCO3⁻ and CO3²⁻ to their initial values)
dydt(39) = 0;  % Fix [HCO3⁻] to its initial value to maintain alkalinity
dydt(40) = 0;  % Fix [H2CO3] to its initial value to maintain alkalinity
dydt(41) = 0;  % Fix [CO3²⁻] to its initial value to maintain alkalinity

        %% Clamp phosphates 
dydt(50) = 0;  % Fix HPO4²⁻ to its initial value to maintain alkalinity
dydt(51) = 0;  % Fix H2PO4⁻ to its initial value to maintain alkalinity


dydt = t_f.(app.t_unit) * dydt;
end

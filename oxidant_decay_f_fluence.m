%% run_oxidant_fluence.m
% This script loops over oxidant conditions (Cl2, NH2Cl, NHCl2, H2O2),
% running the ODE simulation for each and then plotting the normalized concentration
% (C/C₀) versus fluence (mJ/cm^2).

close all; clear; clc;

%% Set basic simulation parameters and water matrix properties
app.t_unit = 'min';             % time unit used in the simulation
% Conversion factors for time:
t_f.sec = 1;                    
t_f.min = 60;  % when time is in minutes, multiply derivatives by 60

% UV irradiance (mW/cm2)
app.E0_mWcm2 = 10;    % incident light intensity in mW/cm2
app.E0_Wm2 = app.E0_mWcm2 * 1; % Energy in W/m2 (J/m2.s)
app.E0 = app.E0_Wm2 * (254e-9 / (6.626e-34 * 3e8 * 6.022e23)); % variable E0 in 'moles of photons/m²·s' | 254nm in m | planck constant in J.s | C in m/S | Avogadro in mol-1

% (Note: 1 mW = 1 mJ/s so time [s] = fluence (mJ/cm2) / (E0_mWcm2))

% pH and temperature
app.NH3_mgL    = 0;       % NH3 (mg N/L)
app.pH = 5.5;                         % pH
app.T_Celcius = 20;                 % temperature in °C
app.T = app.T_Celcius + 273.15;       % in Kelvin

% Water matrix and background properties:
app.alkalinity = 0;  % mg/L as CaCO3
app.DOC        = 0;    % mg C/L
app.chloride   = 0;    % mg/L
app.phosphate  = 2;       % Inorganic Phosphate (mMol/L)
app.sc         = 3e3;  % background scavenging rate (s^-1) [not used here]

% UV source geometry:
app.L = 3.4;         % optical path length (cm)
app.V = 0.0528;    % illuminated volume (L)
app.A = 0.002642;  % illuminated area (m^2)

% Background absorbance:
app.transmittance = 1;
app.a_background = -log10(app.transmittance);

% Other constituents (for example dioxane and carbamazepine)
app.dioxane        = 1e-6;
app.carbamazepine  = 1e-6;
app.DBCP    = 1.83e-6;    % Initial carbamazepine (mol/L)

%% Oxidant conditions:
% The four oxidants are specified via these fields:
%   app.Cl2_mgL   (chlorine, which is partitioned into HOCl and OCl-),
%   app.NH2Cl_mgL,
%   app.NHCl2_mgL,
%   app.H2O2_mgL.
%
% In each simulation, we set one oxidant = 5 mg/L and the others = 0.
% (Conversion factors: see below.)

% Initialize all oxidant fields to zero.
app.Cl2_mgL   = 0;
app.NH2Cl_mgL = 0;
app.NHCl2_mgL = 0;
app.H2O2_mgL  = 0;

oxidants = {'Cl2','NH2Cl','NHCl2','H2O2'};
nOx = length(oxidants);

%% Define fluence axis (mJ/cm^2) and corresponding irradiation time (in minutes)
fluence_array = linspace(0,4500,500);  % 200 points from 0 to 1000 mJ/cm^2
% Time (in seconds) = fluence (mJ/cm^2) / irradiance (mJ/cm^2.s).
% Convert seconds to minutes.
time_array_min = (fluence_array / app.E0_mWcm2) / 60;

%% Set ODE options (as in your run.m)
options = odeset('RelTol',1e-8, 'AbsTol',1e-15, 'NonNegative', 1:63, 'MaxStep', 1e-2);

%% Preallocate matrix to store normalized concentration curves
% Each row corresponds to one oxidant condition.
Cnorm = zeros(nOx, length(fluence_array));

%% Loop over each oxidant scenario
for iOx = 1:nOx
    % Reset all oxidants to 0:
    app.Cl2_mgL   = 0;
    app.NH2Cl_mgL = 0;
    app.NHCl2_mgL = 0;
    app.H2O2_mgL  = 0;
    
    % Set the current oxidant to 5 mg/L:
    switch oxidants{iOx}
        case 'Cl2'
            app.Cl2_mgL = 2.836;
        case 'NH2Cl'
            app.NH2Cl_mgL = 2.836;
        case 'NHCl2'
            app.NHCl2_mgL = 2.836;
        case 'H2O2'
            app.H2O2_mgL = 2.836;
    end
    
    %% Build the initial state vector y0 (based on run.m)
    % Note: Only the oxidant-related species are explicitly shown.
    y0 = zeros(63,1);
    y0(10) = app.NH3_mgL / 14000;         % NH4⁺ (Ammonium ion)
    % NH2Cl (species 1)
    y0(1)  = app.NH2Cl_mgL / 71000;
    % NHCl2 (species 2)
    y0(2)  = app.NHCl2_mgL / (2 * 71000);
    % Cl2 (partitioned into HOCl and OCl-: species 3 and 4)
    y0(3)  = app.Cl2_mgL/(1+10^(app.pH-7.52))/(71*1000);
    y0(4)  = app.Cl2_mgL/(71*1000) - y0(3);
    % H2O2 (species 11)
    y0(11) = app.H2O2_mgL / 34000;
    
    % Other species can be computed as in run.m; here we include a few key ones.
    % (TIC equilibria for carbonate system:)
    H = 10^(-app.pH);
    TIC_molesL = app.alkalinity * 1.22 / 61000;
    K1 = 10^-6.77;
    K2 = 10^-9.93;


    Den = H^2 + H*K1 + K1*K2;
    H2CO3 = TIC_molesL * H^2 / Den;
    HCO3   = TIC_molesL * H*K1 / Den;
    CO3    = TIC_molesL * K1*K2 / Den;
    y0(39) = HCO3;
    y0(40) = H2CO3;
    y0(41) = CO3;

%% Phosphates equilibria at time 0

Phosphate_moles = app.phosphate/1000;

Ka1 = 10^(-2.1);
Ka2 = 10^(-7.2);
Ka3 = 10^(-12.3);

H = 10^(-app.pH);
D = 1 + (Ka1/H) + (Ka1*Ka2)/(H^2) + (Ka1*Ka2*Ka3)/(H^3);

H3PO4  = Phosphate_moles / D;
H2PO4  = Phosphate_moles * (Ka1/H) / D;
HPO4   = Phosphate_moles * (Ka1*Ka2)/(H^2) / D;
PO4    = Phosphate_moles * (Ka1*Ka2*Ka3)/(H^3) / D;

    % pH-dependent species:
    y0(37) = 10^(-app.pH);
    y0(38) = 10^(-14) / y0(37);
    
    % (For dioxane and carbamazepine, we use the values from app.)
    y0(34) = app.dioxane;
    y0(61) = app.carbamazepine;
    y0(62) = app.DBCP;
    % (Other species are set to 0 or computed within merged_ODE.)
    
    %% Determine simulation time span.
    % Run long enough to cover the maximum time corresponding to the highest fluence.
    T_max = max(time_array_min) * 1.2;  % 20% extra margin
    tspan = 0:0.01:T_max;
    
    %% Solve the ODE system using merged_ODE
    [t,y] = ode15s(@(t,y) merged_ODE(t, y, t_f, app), tspan, y0, options);
    
    %% Retrieve the oxidant concentration versus time and compute C/C0.
    % Choose the species index (or combination) as follows:
    switch oxidants{iOx}
        case 'Cl2'
            % Total “Cl2” (chlorine) is taken as the sum of HOCl (species 3) and OCl⁻ (species 4).
            C_t = y(:,3) + y(:,4);
            % The initial concentration is: app.Cl2_mgL/(71*1000)
            C0 = app.Cl2_mgL / (71*1000);
        case 'NH2Cl'
            C_t = y(:,1);
            C0 = app.NH2Cl_mgL / 71000;
        case 'NHCl2'
            C_t = y(:,2);
            C0 = app.NHCl2_mgL / (2*71000);
        case 'H2O2'
            C_t = y(:,11);
            C0 = app.H2O2_mgL / 34000;
    end
    
    %% Instead of interpolating, find the closest time in the solver output.
    for j = 1:length(time_array_min)
        t_target = time_array_min(j);
        % Find the index where the absolute difference is minimal.
        [~, idx] = min(abs(t - t_target));
        C_current = C_t(idx);
        Cnorm(iOx,j) = C_current / C0;
    end
end

%% Plot all oxidant decay curves vs. fluence on one graph
figure;
hold on;
colors = lines(nOx);
for iOx = 1:nOx
    plot(fluence_array, Cnorm(iOx,:), 'Color', colors(iOx,:), 'LineWidth', 2, ...
         'DisplayName', oxidants{iOx});
end
xlabel('Fluence (mJ/cm^2)');
ylabel('C/C_0');
title('Normalized Oxidant Concentration vs. Fluence');
legend('Location','best');
grid on;
hold off;

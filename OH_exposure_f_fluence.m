%% OH_exposure_vs_fluence.m
% This script loops over oxidant conditions (Cl2, NH2Cl, NHCl2, H2O2),
% runs the ODE simulation for each, and plots the cumulative OH exposure
% (integral of [OH]) versus fluence (mJ/cm^2), with curves color coded by oxidant.

close all; clear; clc;

%% Set basic simulation parameters and water matrix properties
app.t_unit = 'min';             % time unit used in the simulation
t_f.sec = 1;                    
t_f.min = 60;  % conversion factor for minutes

% UV irradiance (mW/cm2)
app.E0_mWcm2 = 10;    % incident light intensity in mW/cm2
app.E0_Wm2 = app.E0_mWcm2 * 1; % in W/m2
% Convert energy to mol photons/m^2·s (using 254 nm)
app.E0 = app.E0_Wm2 * (254e-9 / (6.626e-34 * 3e8 * 6.022e23));

% pH and temperature
app.NH3_mgL    = 0;       
app.pH         = 5.6;                         
app.T_Celcius = 20;                 
app.T = app.T_Celcius + 273.15;       

% Water matrix and background properties:
app.alkalinity = 0;  % mg/L as CaCO3
app.DOC        = 0;    
app.chloride   = 0;    
app.phosphate   = 2; %mmol/L
app.sc         = 3e3;  % background scavenging rate (s^-1)

% UV source geometry:
app.L = 3.4;         % optical path length (cm)
app.V = 0.0528;    % illuminated volume (L)
app.A = 0.002642;  % illuminated area (m^2)

% Background absorbance:
app.transmittance = 1;
app.a_background = -log10(app.transmittance);

% Other constituents:
app.dioxane       = 1e-6;
app.carbamazepine = 1e-6;
app.DBCP          = 1.83e-6;

%% Oxidant conditions:
% Initialize oxidant fields to zero.
app.Cl2_mgL   = 0;
app.NH2Cl_mgL = 0;
app.NHCl2_mgL = 0;
app.H2O2_mgL  = 0;

oxidants = {'Cl2','NH2Cl','NHCl2','H2O2'};
nOx = length(oxidants);

%% Define fluence axis (mJ/cm^2) and corresponding irradiation time (in minutes)
fluence_array = linspace(0,2000,1000);  
% time (min) = fluence (mJ/cm^2) / irradiance (mJ/cm^2.s) / 60
time_array_min = (fluence_array / app.E0_mWcm2) / 60;

%% Set ODE options
options = odeset('RelTol',1e-8, 'AbsTol',1e-15, 'NonNegative', 1:63, 'MaxStep', 1e-2);

%% Preallocate matrix to store cumulative OH exposure curves
OH_exp = zeros(nOx, length(time_array_min));

%% Loop over each oxidant scenario
for iOx = 1:nOx
    % Reset oxidant fields:
    app.Cl2_mgL   = 0;
    app.NH2Cl_mgL = 0;
    app.NHCl2_mgL = 0;
    app.H2O2_mgL  = 0;
    
    % Set current oxidant to 5 mg/L:
    switch oxidants{iOx}
        case 'Cl2'
            app.Cl2_mgL = 5;
        case 'NH2Cl'
            app.NH2Cl_mgL = 5;
        case 'NHCl2'
            app.NHCl2_mgL = 5;
        case 'H2O2'
            app.H2O2_mgL = 5;
    end
    
    %% Build the initial state vector (y0)
    % The following mimics your run.m file:
    y0 = zeros(63,1);
    y0(10) = app.NH3_mgL / 14000;         % NH4⁺
    y0(1)  = app.NH2Cl_mgL / 71000;         % NH2Cl
    y0(2)  = app.NHCl2_mgL / (2*71000);      % NHCl2
    % Partition Cl2 into HOCl (species 3) and OCl⁻ (species 4)
    y0(3)  = app.Cl2_mgL/(1+10^(app.pH-7.52))/(71*1000);
    y0(4)  = app.Cl2_mgL/(71*1000) - y0(3);
    y0(11) = app.H2O2_mgL / 34000;          % H2O2
    
    % Carbonate system (from run.m):
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
    
    % Other constituents:
    y0(34) = app.dioxane;
    y0(61) = app.carbamazepine;
    y0(62) = app.DBCP;
    
    %% Determine simulation time span
    T_max = max(time_array_min) * 1.2;  % add a 20% margin
    tspan = 0:0.01:T_max;
    
    %% Solve the ODE system using merged_ODE
    [t,y] = ode15s(@(t,y) merged_ODE(t, y, t_f, app), tspan, y0, options);
    
    %% Compute cumulative OH exposure
    % Extract OH radical concentration (species 12)
    OH_conc = y(:,12);
    % Compute cumulative OH exposure via cumulative trapezoidal integration
    % (units: mol·min/L)
    cumOH = cumtrapz(t, OH_conc);
    
    % Sample cumulative exposure at times in time_array_min
    OH_exp_i = zeros(1, length(time_array_min));
    for j = 1:length(time_array_min)
        t_target = time_array_min(j);
        [~, idx] = min(abs(t - t_target));
        OH_exp_i(j) = cumOH(idx);
    end
    OH_exp(iOx,:) = OH_exp_i;
end

%% Plot cumulative OH exposure vs. fluence for each oxidant
figure;
hold on;
colors = lines(nOx);
for iOx = 1:nOx
    plot(fluence_array, OH_exp(iOx,:), 'Color', colors(iOx,:), 'LineWidth', 2, ...
         'DisplayName', oxidants{iOx});
end
xlabel('Fluence (mJ/cm^2)');
ylabel('Cumulative OH Exposure (mol·min/L)');
title('Cumulative OH Exposure vs. Fluence by Oxidant');
legend('Location','best');
grid on;
hold off;

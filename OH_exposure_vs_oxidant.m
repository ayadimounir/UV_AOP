%% OH_exposure_vs_oxidant.m
% This script computes the cumulative OH exposure (integral of [OH] over time)
% at a fixed fluence of 1000 mJ/cm^2 as a function of the initial oxidant concentration.
% Four oxidants are considered (Cl2, NH2Cl, NHCl2, and H2O2). For each oxidant,
% its initial concentration (in mg/L) is varied over a specified range. The simulation
% uses an irradiance of 0.1 mW/cm^2. For a fluence of 1000 mJ/cm^2, the required exposure time is:
%
%    time (s) = fluence (mJ/cm^2) / irradiance (mW/cm^2)
%             = 1000 / 0.1 = 10,000 s ≈ 166.67 minutes.
%
% The script runs the ODE simulation (using merged_ODE) for each parameter value and then
% plots the resulting cumulative OH exposure versus the initial oxidant concentration.

close all; clear; clc;

%% Simulation Parameters
app.t_unit = 'sec';   % time unit is minutes
t_f.min = 60;         % conversion factor: 1 min = 60 s
t_f.sec = 1;

% UV Source parameters:
app.E0_mWcm2 = 10; % irradiance in mW/cm^2
app.E0_Wm2 = app.E0_mWcm2 * 1; % convert to W/m^2
% Convert to moles of photons/m^2·s (using 254 nm light):
app.E0 = app.E0_Wm2 * (254e-9 / (6.626e-34 * 3e8 * 6.022e23));

% Set desired fluence (mJ/cm^2)
fluence = 1200; 

% Calculate required exposure time (in seconds) and convert to minutes:
T_required_sec = fluence / app.E0_mWcm2;  % 1000/0.1 = 10,000 s
T_required_min = T_required_sec / 60;       % ≈ 166.67 min

% Display the irradiance and exposure time
fprintf('Irradiance: %.2f mW/cm^2\n', app.E0_mWcm2);
fprintf('Exposure time for %.0f mJ/cm^2 fluence: %.2f minutes\n', fluence, T_required_min);

% Water matrix properties:
app.NH3_mgL    = 0;
app.pH = 5.5;
app.T_Celcius = 20;
app.T = app.T_Celcius + 273.15;

app.alkalinity = 0;  % mg/L as CaCO3
app.DOC        = 0;
app.chloride   = 0;
app.phosphate  = 2;       % Inorganic Phosphate (mMol/L)
app.sc         = 3e3;  % background scavenging rate (s^-1)

% UV geometry:
app.L = 3.5;         % optical path length (cm)
app.V = 0.0528;    % illuminated volume (L)
app.A = 0.002642;  % illuminated area (m^2)

% Background absorbance:
app.transmittance = 1;
app.a_background = -log10(app.transmittance);

% Other constituents:
app.dioxane       = 0e-6;
app.carbamazepine = 0e-6;
app.DBCP          = 0e-6;

%% Define oxidant types and concentration range (mg/L)
oxidants = {'Cl2', 'NH2Cl', 'NHCl2', 'H2O2'};
nOx = length(oxidants);
conc_range = linspace(0.5, 10, 30);  % test concentrations from 0 to 10 mg/L (20 points)

% Preallocate matrix to store cumulative OH exposure results 
% (rows: oxidant type, columns: initial concentration values)
OH_exposure = zeros(nOx, length(conc_range));

%% Set ODE options
options = odeset('RelTol',1e-8, 'AbsTol',1e-15, 'NonNegative', 1:63, 'MaxStep', 1e-2);

% Define simulation time span: simulate up to 20% more than T_required_min
T_max = T_required_sec * 1.2;
tspan = 0:0.01:T_max;

%% Loop over each oxidant type and concentration value
for iOx = 1:nOx
    for iConc = 1:length(conc_range)
        % Reset all oxidant fields:
        app.Cl2_mgL   = 0;
        app.NH2Cl_mgL = 0;
        app.NHCl2_mgL = 0;
        app.H2O2_mgL  = 0;
        
        % Set the current oxidant to the test concentration
        testConc = conc_range(iConc);
        switch oxidants{iOx}
            case 'Cl2'
                app.Cl2_mgL = testConc;
            case 'NH2Cl'
                app.NH2Cl_mgL = testConc;
            case 'NHCl2'
                app.NHCl2_mgL = testConc;
            case 'H2O2'
                app.H2O2_mgL = testConc;
        end
        
        %% Build the initial state vector y0 (following run.m conventions)
        y0 = zeros(63,1);
        % Ammonium (species 10)
        y0(10) = app.NH3_mgL / 14000;
        % Chloramine species:
        y0(1)  = app.NH2Cl_mgL / 71000;        % NH2Cl (species 1)
        y0(2)  = app.NHCl2_mgL / (2*71000);       % NHCl2 (species 2)
        % For Cl2, partition into HOCl (species 3) and OCl⁻ (species 4)
        y0(3)  = app.Cl2_mgL/(1+10^(app.pH-7.52))/(71*1000);
        y0(4)  = app.Cl2_mgL/(71*1000) - y0(3);
        % H2O2 (species 11)
        y0(11) = app.H2O2_mgL / 34000;
        
        % Carbonate system equilibria:
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
        
        % Other species:
        y0(34) = app.dioxane;
        y0(61) = app.carbamazepine;
        y0(62) = app.DBCP;
        
        %% Solve the ODE system (using merged_ODE)
        [t, y] = ode15s(@(t,y) merged_ODE(t, y, t_f, app), tspan, y0, options);
        
        %% Compute cumulative OH exposure
           %% Extract cumulative OH exposure directly from y(63)
        OH_val = interp1(t, y(:,63), T_required_sec, 'linear', 'extrap');
        OH_exposure(iOx, iConc) = OH_val;
    end
end

%% Plot the results
figure;
hold on;
markers = {'-o','-s','-^','-d'};
for iOx = 1:nOx
    plot(conc_range, OH_exposure(iOx,:), markers{iOx}, 'LineWidth', 2, ...
         'DisplayName', oxidants{iOx});
end
xlabel('Initial Oxidant Concentration (mg/L)');
ylabel('Cumulative OH Exposure (from y(63))');
title(sprintf('OH Exposure at %.0f mJ/cm^2 Fluence\n(E0 = %.2f mW/cm^2, t = %.2f min)', ...
      fluence, app.E0_mWcm2, T_required_min));
legend('Location','best');
grid on;
hold off;

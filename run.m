% run.m
close all;
clear;
clc;

%% Define unit, temperature, lamp intensity, etc.
app.t_unit = 'min';               % Choose time unit: 'sec' or 'min'
tspan = 0:0.01:25;                 % Time range in chosen unit
app.t_switch = 20;                 % UV light turns on at t = 2 minutes

app.T_Celcius = 20;               % Temperature in °C (will be converted to Kelvin)

% Time conversion factors:
t_f.sec = 1;
t_f.min = 60;

%% User-defined oxidants (mg/L)
app.Cl2_mgL   = 0;    % ClO- (mg Cl2/L)
app.NH2Cl_mgL = 0;   % NH2Cl (mg Cl2/L)
app.NHCl2_mgL = 4;    % NHCl2 (mg Cl2/L)
app.H2O2_mgL  = 0;    % H2O2 (mg/L)
app.NCl3_mgL  = 0;    % NCl3 (mg Cl2/L)

%% Water matrix initial parameters
app.NH3_mgL    = 0;       % NH3 (mg N/L)
app.pH         = 7;      % pH
app.sc         = 5e5;     % Background scavenging (s^-1) 

%% UV Source
% % We either use energy based (E0_mWm2) or photon based (E0)

app.E0_mWcm2 = 10; % Energy in mW/cm2 (mJ/cm2.s)
app.E0_Wm2 = app.E0_mWcm2 * 1; % Energy in W/m2 (J/m2.s)
app.E0 = app.E0_Wm2 * (254e-9 / (6.626e-34 * 3e8 * 6.022e23)); % variable E0 in 'moles of photons/m²·s' | 254nm in m | planck constant in J.s | C in m/S | Avogadro in mol-1


app.L = 3.5;        % Fixed optical path length (cm)
app.V = 0.0528;   % Illuminated volume (L)
app.A = 0.002642;   % Illuminated Area (m^2)

%% water matrix parameters

app.alkalinity = 45;      % Alkalinity (mg/L as CaCO3)
app.DOC        = 0;       % DOC (mg C/L)
app.chloride   = 40;      % Cl- (mg/L)
app.phosphate  = 0;       % Inorganic Phosphate (mMol/L)

app.dioxane    = 1e-6;    % Initial 1,4-Dioxane (mol/L)
app.caffeine    = 1e-6;    % Initial caffeine (mol/L)
app.sucralose    = 1e-6;    % Initial sucralose (mol/L)
app.carbamazepine    = 1e-6;    % Initial carbamazepine (mol/L)


%background absorbance (excluding oxidants)
app.transmittance = 0.717; % Example transmittance value (99%)

% Calculate Absorbance (A) from Transmittance (T)
app.a_background = -log10(app.transmittance);


% Convert temperature from °C to Kelvin
app.T = app.T_Celcius + 273.15;

%% Initialize species concentrations (state vector of 62×1, in mol/L)

%% TIC equilibria at time 0

H=10^(-app.pH); %% variables here are used like in anonymous functions, they shouldn't be used somewhere else unless good reason
TIC_molesL = app.alkalinity*1.22/61000;

K1 = 10^-6.77; % H2CO3 <-> HCO3 + H
    K2 = 10^-9.93; % HCO3 <-> CO3 + H
    Den = H*H + H*K1 + K1*K2;
    H2CO3 = TIC_molesL*H*H/Den;
    HCO3 = TIC_molesL*H*K1/Den;
    CO3 = TIC_molesL*K1*K2/Den;

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
%% Initial state vector build-up

y0 = zeros(64,1);

% Chloramine species
y0(1)  = app.NH2Cl_mgL / 71000;        % NH2Cl (monochloramine)
y0(2)  = app.NHCl2_mgL / (2 * 71000);  % NHCl2 (dichloramine)
y0(3)  = app.Cl2_mgL/(1+10^(app.pH-7.52))/(71*1000); % HOCl (hypochlorous acid)
y0(4)  = app.Cl2_mgL/(71 * 1000) - y0(3);  % OCl⁻ (hypochlorite ion)




% Nitrogen compounds
y0(8)  = app.NCl3_mgL / (3 * 71000);  % NCl3 (trichloramine)
y0(9)  = 0;                           % NH3 (Ammonia)
y0(10) = app.NH3_mgL / 14000;         % NH4⁺ (Ammonium ion)

% H2O2
y0(11)  = app.H2O2_mgL/(34 * 1000);  % OCl⁻ (hypochlorite ion)


% Organic pollutants
y0(34) = app.dioxane;  % 1,4-Dioxane
y0(35) = 0;            % DOC (Dissolved organic carbon)

% Oxygen species
y0(36) = 0;  % O2 (Molecular oxygen)

% pH-dependent species
y0(37) = 10^(-app.pH);        % H⁺ (Hydrogen ion concentration from pH)
y0(38) = 10^(-14)/y0(37);     % OH⁻ (Hydroxide ion from water equilibrium)

% Carbonate species
y0(39) = HCO3;   % HCO3⁻ (Bicarbonate ion)
y0(40) = H2CO3;  % H2CO3 (Carbonic acid)
y0(41) = CO3;    % CO₃²⁻ (Carbonate ion)

% Chloride-related species
y0(42) = app.chloride / (35.5 * 1000);  % Cl⁻ (Chloride ion)

% phosphate species

y0(51) = H2PO4;  % H2PO4-
y0(50) = HPO4;   % HPO4^2-

% Pharmaceuticals and emerging contaminants
y0(61) = app.caffeine;  % caffeine
y0(62) = app.sucralose;  % sucralose
y0(63) = 0;
y0(64) = app.carbamazepine; % Carbamazepine


%% Set ODE options: update NonNegative indices to include species 1:60
options = odeset('RelTol',1e-8, 'AbsTol',1e-15, 'NonNegative', 1:64, 'MaxStep', 1e-2);

%% Solve the ODE system:
[t, y] = ode15s(@(t,y) merged_ODE(t, y, t_f, app), tspan, y0, options);

%% Determine time necessary for 0.5-log elimination
% Calculate the 0.5-log thresholds for 1,4-dioxane, caffeine, and sucralose
threshold_diox = y0(34) * 10^(-0.5);
threshold_carb = y0(61) * 10^(-0.5);
threshold_sucralose = y0(62) * 10^(-0.5);

% Initialize elimination times as NaN (in case threshold is not reached)
time_diox = NaN;
time_carb = NaN;
time_sucralose = NaN;

% Find the crossing time for 1,4-dioxane (species 34)
for i = 1:length(t)-1
    if (y(i,34) > threshold_diox && y(i+1,34) <= threshold_diox) || ...
       (y(i,34) < threshold_diox && y(i+1,34) >= threshold_diox)
        time_diox = interp1([y(i,34), y(i+1,34)], [t(i), t(i+1)], threshold_diox);
        break;
    end
end

% Find the crossing time for caffeine (species 61)
for i = 1:length(t)-1
    if (y(i,61) > threshold_carb && y(i+1,61) <= threshold_carb) || ...
       (y(i,61) < threshold_carb && y(i+1,61) >= threshold_carb)
        time_carb = interp1([y(i,61), y(i+1,61)], [t(i), t(i+1)], threshold_carb);
        break;
    end
end

% Find the crossing time for sucralose (species 62)
for i = 1:length(t)-1
    if (y(i,62) > threshold_sucralose && y(i+1,62) <= threshold_sucralose) || ...
       (y(i,62) < threshold_sucralose && y(i+1,62) >= threshold_sucralose)
        time_sucralose = interp1([y(i,62), y(i+1,62)], [t(i), t(i+1)], threshold_sucralose);
        break;
    end
end

% Create the output object as a struct:
time_for_0_5log_elimination = struct('dioxane', time_diox, ...
                                     'caffeine', time_carb, ...
                                     'sucralose', time_sucralose);


%% Plot the results using a separate function
    plot_results(t, y, app);

%% extract OH exposure
OHexposure = y(:,64); % Extract OH exposure values at each time step
function plot_results(t, y, app)
%% Plot results using tiledlayout for improved vertical space

% --- Global figure/axes defaults ---
set(groot, 'defaultAxesXGrid','on', ...
           'defaultAxesYGrid','on', ...
           'defaultAxesBox','on', ...
           'defaultAxesFontName','Arial', ...
           'defaultAxesFontSize',12, ...
           'defaultAxesLabelFontSizeMultiplier',1.2, ...
           'defaultAxesTitleFontSizeMultiplier',1.2);



figure;
% Create a 3x3 tiled layout with compact spacing
tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact');


% --- Read experimental data ---
expData = readtable('Rayox_complete_dataset.xlsx', ...
    'Sheet', 'sheet1', 'VariableNamingRule','preserve');

% --- Apply filters ---
expData = expData( ...
    expData.('pHi') == 7 & ...
    expData.('NH4 C0') > 0 & ...
    strcmp(expData.('Type'), 'Wastewater') & ...
    expData.('Cl2 injected (mgCl2/L)') == 4, :);

% --- Extract columns after filtering ---
t_exp         = expData.('total time(minute)');
freeCl_exp    = expData.('Cl2free (mgCl2/L)');
totalCl_exp   = expData.('Cl2total (mgCl2/L)'); 
dioxane_exp   = expData.('dioxane C/C0');
sucralose_exp = expData.('sucralose C/C0');
caffeine_exp  = expData.('caffeine C/C0');
carbamazepine_exp  = expData.('Carbamazepine C/C0');
pH_exp        = expData.('pH');


% --- Extract condition info for subplot titles ---
pHi_val  = expData.('pHi')(1);
Cl2_val  = expData.('Cl2 injected (mgCl2/L)')(1);
title_suffix = sprintf('(wastewater, pH = %.1f, Cl_2 = %.1f mg/L)', pHi_val, Cl2_val);


% Identify unique pH values and assign a color to each
unique_pH = unique(pH_exp);
colors = lines(length(unique_pH));  % Or manually define colors if preferred


%-------------------------------
% Subplot 1: Oxidants Concentration with Experimental Free Chlorine Data
%-------------------------------
nexttile;
grid on;
hold on;

% Molar masses in mg/mol (from run.m)
MM_Cl2 = 71000; % For all chlorine species (mg Cl2/mol)
MM_H2O2 = 34000; % For H2O2 (mg/mol)

% Model species (converted from mol/L to mg/L for plotting)
plot(t, y(:,3) * MM_Cl2, 'k', ...      % [HOCl] in mg/L as Cl2
     t, y(:,4) * MM_Cl2, 'r', ...      % [OCl-] in mg/L as Cl2
     t, y(:,1) * MM_Cl2, 'g', ...      % [NH2Cl] in mg/L as Cl2
     t, y(:,2) * 2 * MM_Cl2, 'b', ...  % [NHCl2] in mg/L as Cl2
     t, y(:,11) * MM_H2O2, 'y', ...    % [H2O2] in mg/L
     'LineWidth', 1.2);

% Model's total Chlorine = HOCl + OCl- + NH2Cl + NHCl2 (all converted to mg/L)
plot(t, ((y(:,3) + y(:,4)) * MM_Cl2 + y(:,1) * MM_Cl2 + y(:,2) * 2 * MM_Cl2), '--', ...
     'Color', [0 1 1], ...
     'LineWidth', 1.5);

% Model's free Chlorine = HOCl + OCl- + NH2Cl + NHCl2 (all converted to mg/L)
plot(t, ((y(:,3) + y(:,4)) * MM_Cl2), '--', ...
     'Color', [1 0 0], ...
     'LineWidth', 1.5);

% Experimental free chlorine (solid red circles)
% Note: Multiplied by initial concentration to convert C/C0 to mg/L
plot(t_exp, freeCl_exp, 'o', ...
     'MarkerEdgeColor', 'r', ...
     'MarkerFaceColor', 'r', ...
     'LineWidth', 1.5, ...
     'DisplayName', 'Exp Free Cl');

% Experimental total chlorine (solid cyan circles)
% Note: Multiplied by initial concentration to convert C/C0 to mg/L
plot(t_exp, totalCl_exp, 'o', ...
     'MarkerEdgeColor', 'c', ...
     'MarkerFaceColor', 'c', ...
     'LineWidth', 1.5, ...
     'DisplayName', 'Exp Total Cl');

hold off;

% Legend & labels
legend( ...
    {'HOCl', 'OCl^-', 'NH2Cl', 'NHCl2', 'H2O2', 'Total Cl', 'Free Cl', 'Exp Free Cl', 'Exp Total Cl'}, ...
    'Location', 'best' );
xlabel(['Time (' app.t_unit ')']);
ylabel('Concentration (mg/L)');
title(['Oxidant Concentrations ' title_suffix]);


% Subplot 2: Ammonia Concentration (mgN/L)
nexttile;
grid on; box on;
% Replace k with the actual column index for ammonia (NH3/NH4+) in y
ammonia_mol = y(:, 10);  

% Convert mol/L to mgN/L
MM_N = 14000;  % mg N per mol N
ammonia_mgN_L = ammonia_mol * MM_N;

plot(t, ammonia_mgN_L, 'b', 'LineWidth', 1.2);
xlabel(['Time (' app.t_unit ')']);
ylabel('mg N / L');
title('Ammonia Concentration ');
legend('NH_3 / NH_4^+','Location','best');
ylim([0 1.1]); 



% Subplot 3: pH Over Time
nexttile;
grid on; box on;
pH = -log10(y(:,37));  % pH computed from H+ concentration
plot(t, pH, 'b');
xlabel(['Time (' app.t_unit ')']);
ylabel('pH');
title('pH Over Time');
ylim([3 11]);  % Set pH axis to range between 3 and 11

% -------------------------------
% Subplot 4: Micropollutant Concentrations
% -------------------------------
nexttile;
grid on; box on;

% 0.5‑log reduction reference value
log_line = 10^(-0.5);

% Plot computed micropollutants
plot(t, y(:,34)/y(1,34), 'b',  'DisplayName','Model 1,4‑Dioxane'); hold on;
plot(t, y(:,61)/y(1,61), 'g',  'DisplayName','Model Caffeine');
plot(t, y(:,62)/y(1,62), 'r',  'DisplayName','Model Sucralose');
plot(t, y(:,64)/y(1,64), 'y',  'DisplayName','Model Carbamazepine');

% 0.5‑log reduction reference line
plot(t, ones(size(t))*log_line, '--k', 'DisplayName','0.5-log reduction');

% Now plot experimental points with solid markers
% Dioxane (blue circles)
plot(t_exp, dioxane_exp, 'o', ...
     'MarkerEdgeColor','b', 'MarkerFaceColor','b', ...
     'DisplayName','Exp Dioxane');
% Caffeine (green squares)
plot(t_exp, caffeine_exp, 's', ...
     'MarkerEdgeColor','g', 'MarkerFaceColor','g', ...
     'DisplayName','Exp Caffeine');
% Sucralose (red triangles)
plot(t_exp, sucralose_exp, '^', ...
     'MarkerEdgeColor','r', 'MarkerFaceColor','r', ...
     'DisplayName','Exp Sucralose');
% carbamazepine (dark yellow)
plot(t_exp, carbamazepine_exp, '^', ...
     'MarkerEdgeColor','y', 'MarkerFaceColor','y', ...
     'DisplayName','Exp carbamazepine');

hold off;
xlabel(['Time (' app.t_unit ')']);
ylabel('C/C_0');
title(['Target micropollutant ' title_suffix]);
legend('show','Location','best');


% Subplot 5: [·OH] Exposure
nexttile;
grid on; box on;
plot(t, y(:,63), 'k');   % [·OH] exposure
legend('*OH exposure','Location','best');
xlabel(['Time (' app.t_unit ')']);
ylabel('mol/L');
title('[·OH] Exposure');

% Subplot 6: Carbonate Species Concentrations
nexttile;
grid on; box on;
plot(t, y(:,40), 'k', ...   % H2CO3
     t, y(:,39), 'r', ...   % HCO3⁻
     t, y(:,41), 'b');      % CO3²⁻
legend('H2CO3','HCO3^-','CO3^{2-}','Location','best');
xlabel(['Time (' app.t_unit ')']);
ylabel('mol/L');
title('Carbonate Species Concentrations');

% Subplot 7: Phosphate Species Concentrations
nexttile;
grid on; box on;
plot(t, y(:,50), 'm', ...   % HPO4²⁻
     t, y(:,51), 'c');      % H2PO4⁻
legend('HPO4^{2-}','H2PO4^-','Location','best');
xlabel(['Time (' app.t_unit ')']);
ylabel('mol/L');
title('Phosphate Species Concentrations');

% Subplot 8: Total UV Transmittance
nexttile;
grid on; box on;
total_absorbance = 62 * y(:,3) + 66 * y(:,4) + 371 * y(:,1) + ...
                   142 * y(:,2) + 18 * y(:,11) + 228 * y(:,13) + app.a_background;
UVT = 10.^(-total_absorbance) * 100;
plot(t, UVT, 'k');
legend('UVT','Location','best');
xlabel(['Time (' app.t_unit ')']);
ylabel('UVT (%)');
title('Total UV Transmittance');

% Subplot 9: Total Absorbance
nexttile;
grid on; box on;
plot(t, total_absorbance, 'k');
legend('Total Absorbance','Location','best');
xlabel(['Time (' app.t_unit ')']);
ylabel('Absorbance');
title('Total Absorbance');

end

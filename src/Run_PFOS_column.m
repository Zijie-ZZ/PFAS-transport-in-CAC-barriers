% run_PFOS_column.m

% Austhors: Zijie Zheng and Marc Hesse
% The University of Texas at Austin

% This script reproduces PFOS column breakthrough curves using a
% 1D advection–dispersion model with nonlinear Freundlich sorption
% and a two-site (equilibrium + kinetic) formulation.
%
% HOW TO REUSE THIS SCRIPT:
% - To simulate a different PFAS (e.g., PFOA, PFHxS), modify the
%   sorption and kinetic parameters in the "Model parameters" section.
% - To simulate a different column or flow condition, modify the
%   injection duration and PV conversion parameters.
% - The core transport physics is implemented in
%   ADE_Freundlich_Kinetic_ana_function.m (do not edit unless needed).
%
% Users are encouraged to treat this script as a driver/wrapper
% around the core solver for their own simulations.

clear; close all; clc
% Always use paths relative to this script location so the code
% runs regardless of where the repository is downloaded
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));   % assumes script is in repoRoot/src
addpath(fullfile(repoRoot,'src'));


%% ---- Load experimental data ----
% This function loads PFOS column breakthrough data and returns a
% clean structure array Exp(i) with fields:
%   - Q_mL_h  : flow rate
%   - t_days : time vector (days)
%   - C_C0   : normalized effluent concentration
%
% To use your own data:
%   - Replace PFOS_data_plot.mat in the data/ folder, OR
%   - Modify load_PFOS_data_from_mat.m to read a different format (CSV, MAT)

dataFile = fullfile(repoRoot,'data','PFOS_data_plot.mat');
Exp = load_PFOS_data_from_mat(dataFile);   % returns Exp(i) with Q, t_days, C_C0

%% ---- Model parameters (PFAS-specific; edit to simulate other compounds) ----
% Influent concentration
par.c0    = 0.20;      % g/m^3

% Freundlich sorption parameters (compound-specific)
par.Kf    = 315.09;    % (g/kg)/(g/m^3)^nf
par.nf    = 0.835;     % [-]

% Bulk density of sorbent/amended media
par.rho_b = 0.0157;    % (unit must match your solver)

% Kinetic mass-transfer rate constant
par.alpha = 1.046;     % 1/day

% Fraction of sorption sites at equilibrium:
%   f < 1  : mixed equilibrium + kinetic sorption
%   f = 1  : No kinetic sorption (instantaneous equilibrium)
f_list = [0.176, 1.0];

% Injection duration (pulse length) for each flow rate. 
% {12,24,36} units (ml/h) - {0.1111, 0.0556, 0.0370} units (m/day)
% These values reflect the experimental injection protocol.
%
% To simulate a different pulse length or continuous injection,
% modify these values accordingly.
tInjectDays_map = containers.Map({12,24,36}, {0.1111, 0.0556, 0.0370}); 


% Total simulation time (days). - adjust if needed
% Should be long enough to capture full breakthrough and elution.
% Here we use 1 day to ensure all cases reach ~50 PV.
tMaxDays = 1.0;

%% ---- PV conversion ----
% pvph_map defines pore volumes per hour for each flow rate,
% computed from column geometry and porosity.
%
% Users may:
% - Modify pvph values to match a different column
% - Set usePV = false to plot results in time instead of PV
pvph_map = containers.Map({12,24,36}, {2.94, 7.33, 10.9});  % PV per hour
usePV = true;

%% ---- Plot ----
figure('Position',[60 120 1350 420]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for iQ = 1:numel(Exp)
    Q = Exp(iQ).Q_mL_h;

    nexttile; hold on; box on

    % Experiment
    xExp = Exp(iQ).t_days;
    if usePV
        xExp = time_days_to_PV(xExp, pvph_map(Q));
    end
    plot(xExp, Exp(iQ).C_C0, 'bo', 'MarkerFaceColor','none', 'MarkerSize',7, ...
        'DisplayName','experiment');

    % Model curves
    for jf = 1:numel(f_list)
        par.f = f_list(jf);
        tInjectDays = tInjectDays_map(Q);
       
        % Call the core transport solver
        % (users typically do NOT need to modify this function)
        [t_days, Cout] = simulate_single_PFOS_column(par, Q, tInjectDays, tMaxDays);

        xMod = t_days;
        if usePV
            xMod = time_days_to_PV(xMod, pvph_map(Q));
        end

        if par.f < 1
            plot(xMod, Cout/par.c0, 'b-', 'LineWidth',2.2, 'DisplayName','model (best-fit f)');
        else
            plot(xMod, Cout/par.c0, 'r--', 'LineWidth',2.2, 'DisplayName','model (f=1)');
        end
    end

    if usePV
        xlabel('Pore volumes (PV)');
        xlim([0 50]);
    else
        xlabel('Time (days)');
    end

    ylabel('C/C_0');
    title(sprintf('PFOS | Q = %d mL/h', Q));
    ylim([0 1.2]);
    legend('Location','southeast','Box','off');
    set(gca,'FontSize',12);
end

% End of script.
% This driver is intentionally kept simple so it can be adapted
% to other PFAS compounds, column geometries, and flow conditions
% with minimal modification.

%% ===== master driver =====
clear; clc; close all;

project_root = '/Users/luminjhe/Downloads/CHsolvers_package';
addpath(fullfile(project_root,'CahnHilliard_MATLAB_solvers'));
addpath(fullfile(project_root,'src'));

%% ---------- IC ----------
IC_dir = fullfile(project_root,'IC');
% icFile = fullfile(IC_dir,'initial_phi_256_smooth_n_relax_4_from512.csv');
icFile = fullfile(IC_dir,'circle_R01_N256.csv');
phi0   = readmatrix(icFile);

%% ---------- global settings ----------
Tfinal   = 0.001;        % ←  mandatory final time
boundary = 'neumann';
dtList   = logspace(-3,-6,4);     % 1e-3 … 1e-6

% %% ---------- η scan ----------
% etaList = [0, 0.95, 1];
% gamma0_fixed = 4;
% result_eta = run_eta_scan(phi0,dtList,etaList,...
%                           Tfinal,boundary,gamma0_fixed);

%% ---------- γ₀ scan ----------
gammaList = [0,1,2,4];
eta_fixed = 0.95;
result_gamma = run_gamma_scan(phi0,dtList,gammaList,...
                              Tfinal,boundary,eta_fixed);

%% ---------- plotting ----------
% util_plot(result_eta,'eta');
util_plot(result_gamma,'gamma0');

save(fullfile(project_root,'SAV_scan_results.mat'), ...
     'result_eta','result_gamma','dtList', ...
     'etaList','gammaList','Tfinal','-v7.3');
disp('>>>  finished – data written to SAV_scan_results.mat');
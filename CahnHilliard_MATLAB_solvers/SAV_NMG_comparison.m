%% ------------------------------------------------------------
%  CH: NMG  vs  SAV  –  Convergence & Efficiency   (rev-A)
%  *   slope fitted only on non-zero errors
%  *   extra finest step 5e-6 to serve as reference
%  *   initial condition loaded from CSV (128×128)
%  ------------------------------------------------------------
clear;  clc;

%% 1. Load initial condition
file_path = '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_128_smooth_n_relax_4.csv';
phi0      = readmatrix(file_path);          % size 128×128 (double)

%% 2. User parameters
tFinal   = 2e-2;                            % target time 0.02
dtList   = [1e-3 5e-4 2e-4 1e-4 5e-5 2.5e-5 1e-5 5e-6];  % ← added 5e-6
boundary = "neumann";
savePlots = true;

gamma0=4; eta=0.95; xi_flag=1; C0=0;       % SAV parameters

%% 3. Pre-allocate
nCases = numel(dtList);
phiNMG = cell(1,nCases);  phiSAV = cell(1,nCases);
errNMG = zeros(1,nCases); errSAV = zeros(1,nCases);
cpuNMG = zeros(1,nCases); cpuSAV = zeros(1,nCases);

%% 4. Main loop – single output at tFinal
for k = 1:nCases
    dt      = dtList(k);
    nSteps  = ceil(tFinal/dt);
    dt_out  = nSteps;                       % output only at final time

    % ---------- NMG ----------
    tic;
    [~, phi_t, ~, ~] = CahnHilliard_NMG( ...
        phi0, 't_iter', nSteps, 'dt', dt, ...
        'dt_out', dt_out, 'boundary', boundary);
    cpuNMG(k) = toc;
    phiNMG{k} = phi_t(:,:,end);

    % ---------- SAV ----------
    tic;
    [~, phi_t, ~, ~] = CahnHilliard_SAV( ...
        phi0, 't_iter', nSteps, 'dt', dt, ...
        'dt_out', dt_out, 'boundary', boundary, ...
        'gamma0',gamma0,'eta',eta,'xi_flag',xi_flag,'C0',C0);
    cpuSAV(k) = toc;
    phiSAV{k} = phi_t(:,:,end);

    fprintf("Δt = %-9g | NMG %.2fs | SAV %.2fs\n", dt, cpuNMG(k), cpuSAV(k));
end

%% 5. Error vs reference (finest dt = 5e-6)
phiRef_NMG = phiNMG{end};
phiRef_SAV = phiSAV{end};
for k = 1:nCases
    errNMG(k) = norm(phiNMG{k}(:)-phiRef_NMG(:),2);
    errSAV(k) = norm(phiSAV{k}(:)-phiRef_SAV(:),2);
end

%% 6. Remove zero-error point when fitting slopes
mask = errNMG>0 & errSAV>0;   % logical index (last point excluded)

pNMG = polyfit(log(dtList(mask)), log(errNMG(mask)), 1);
pSAV = polyfit(log(dtList(mask)), log(errSAV(mask)), 1);
orderNMG = -pNMG(1);
orderSAV = -pSAV(1);

%% 7. Plot: error vs Δt
figure(1); clf; hold on;
loglog(dtList,errNMG,'-o','LineWidth',1.6);
loglog(dtList,errSAV,'-s','LineWidth',1.6);
set(gca,'XDir','reverse'); grid on;
xlabel('\Delta t'); ylabel('L_2 error (vs.\ finest)');
legend({'NMG','SAV'},'Location','southwest');
title('Temporal convergence at t = 0.02');

txt = sprintf('order:  NMG \\approx %.2f   SAV \\approx %.2f',orderNMG,orderSAV);
text(dtList(4),errNMG(4),txt,'FontSize',9);   % anchor roughly mid-plot

if savePlots, saveas(1,'CH_dt_convergence_loglog_fixed.png'); end

%% 8. Plot: CPU-time vs error
figure(2); clf; hold on;
loglog(cpuNMG,errNMG,'-o','LineWidth',1.6);
loglog(cpuSAV,errSAV,'-s','LineWidth',1.6);
grid on;
xlabel('CPU time (s)'); ylabel('L_2 error');
legend({'NMG','SAV'},'Location','southwest');
title('Efficiency: CPU time vs accuracy');
if savePlots, saveas(2,'CH_cpu_vs_error_fixed.png'); end
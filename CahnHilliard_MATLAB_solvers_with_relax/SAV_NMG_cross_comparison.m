%% ================================================================
%  CH: NMG vs SAV – Cross-convergence, Efficiency, Field Overlay
%  Author : Min-Jhe Lu   (2025-05-08, colormap edition)
% ================================================================

clear;  clc;  close all;

%% 1. Load initial condition ---------------------------------------------
icPath = '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_128_smooth_n_relax_4.csv';
phi0   = readmatrix(icPath);          % 128 × 128 

%% 2. Global parameters ---------------------------------------------------
tFinal   = 2e-2;                                             % target time
dtList   = [1e-3 5e-4 2e-4 1e-4 5e-5 2.5e-5 1e-5 5e-6];      % time-steps
boundary = "neumann";

% SAV parameters
gamma0 = 4;  eta = 0.95;  xi_flag = 1;  C0 = 0;

savePlots = true;
dtTarget  = 1e-4;                                            % visualisation step

%% 3. Storage allocation --------------------------------------------------
n = numel(dtList);
phiNMG = cell(1,n);  phiSAV = cell(1,n);
errNMG = zeros(1,n); errSAV = zeros(1,n); errGap = zeros(1,n);
cpuNMG = zeros(1,n); cpuSAV = zeros(1,n);

%% 4. Run both solvers ----------------------------------------------------
for k = 1:n
    dt     = dtList(k);
    nSteps = ceil(tFinal/dt);
    dt_out = nSteps;                        % single snapshot

    % ----- NMG
    tic;
    [~, out, ~, ~] = CahnHilliard_NMG(phi0, ...
        't_iter',nSteps,'dt',dt,'dt_out',dt_out,'boundary',boundary);
    cpuNMG(k) = toc;
    phiNMG{k} = out(:,:,end);

    % ----- SAV
    tic;
    [~, out, ~, ~] = CahnHilliard_SAV(phi0, ...
        't_iter',nSteps,'dt',dt,'dt_out',dt_out,'boundary',boundary, ...
        'gamma0',gamma0,'eta',eta,'xi_flag',xi_flag,'C0',C0);
    cpuSAV(k) = toc;
    phiSAV{k} = out(:,:,end);

    fprintf('Δt = %-8g | NMG %.2fs | SAV %.2fs\n', dt, cpuNMG(k), cpuSAV(k));
end

%% 5.  Errors relative to finest NMG -------------------------------------
phiRef = phiNMG{end};               % dt = 5e-6 reference

for k = 1:n
    errNMG(k) = norm(phiNMG{k}(:) - phiRef(:), 2);
    errSAV(k) = norm(phiSAV{k}(:) - phiRef(:), 2);
    errGap(k) = norm(phiSAV{k}(:) - phiNMG{k}(:), 2);
end

%% 6.  Convergence order (SAV → ref) -------------------------------------
mask      = errSAV > 0;                          % skip zero point
coeff     = polyfit(log(dtList(mask)), log(errSAV(mask)), 1);
orderSAV  = -coeff(1);

%% 7.  Figure 1 – convergence curves -------------------------------------
figure(1); clf; hold on;
loglog(dtList,errNMG,'-o','LineWidth',1.6);
loglog(dtList,errSAV,'-s','LineWidth',1.6);
loglog(dtList,errGap,'-.^','LineWidth',1.6);
set(gca,'XDir','reverse'); grid on;
xlabel('\Delta t'); ylabel('L_2 error');
legend({'NMG vs ref','SAV vs ref','SAV vs NMG'},'Location','southwest');
title({'Convergence toward finest-step NMG', 't = 0.02'});
text(dtList(4),errSAV(4),sprintf('order (SAV→ref) ≈ %.2f',orderSAV),...
     'FontSize',9);

if savePlots, saveas(1,'CH_cross_convergence.png'); end

%% 8.  Figure 2 – CPU time vs error --------------------------------------
figure(2); clf; hold on;
loglog(cpuNMG,errNMG,'-o','LineWidth',1.6);
loglog(cpuSAV,errSAV,'-s','LineWidth',1.6);
grid on;
xlabel('CPU time (s)'); ylabel('L_2 error vs ref');
legend({'NMG','SAV'},'Location','southwest');
title('Efficiency: CPU time vs accuracy');

if savePlots, saveas(2,'CH_cpu_vs_error.png'); end

%% 9.  Figure 3 – side-by-side fields with red-blue map -------------------
[~,idxCmp] = min(abs(dtList - dtTarget));    % chosen Δt index
phiN = phiNMG{idxCmp};
phiS = phiSAV{idxCmp};
dPhi = phiN - phiS;

% build 1000-level diverging colormap
if exist('redbluecmap','file')
    rbMap = interp1(1:100:1100, redbluecmap, 1:1001);
else
    error('redbluecmap.m not found on MATLAB path.');
end
cLim = [min(phi0(:)) max(phi0(:))];          % shared colour limits

figure('Name','Field comparison','Position',[100 100 1200 400]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% -- NMG
nexttile;
imagesc(phiN); axis equal tight off;
colormap(gca, rbMap); caxis(cLim); colorbar;
title(sprintf('NMG  \\Delta t = %.1e', dtList(idxCmp)));

% -- SAV
nexttile;
imagesc(phiS); axis equal tight off;
colormap(gca, rbMap); caxis(cLim);
title(sprintf('SAV  \\Delta t = %.1e', dtList(idxCmp)));

% -- Difference
nexttile;
imagesc(dPhi); axis equal tight off;
colormap(gca,'jet'); colorbar;
title('\phi_{NMG} - \phi_{SAV}');

if savePlots
    saveas(gcf, sprintf('CH_field_comparison_dt%.0e.png', dtList(idxCmp)));
end
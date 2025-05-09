%% =============================================================
%  Triplet storyboard at ultra-fine steps
%  Δt list: 5e-6  →  2e-7  →  1e-7
% =============================================================

clear; clc; close all;

%% 1.  Initial field ------------------------------------------------------
icPath = '/Users/luminjhe/Downloads/IC/initial_phi_128_smooth_n_relax_4.csv';
phi0   = readmatrix(icPath);

%% 2.  Simulation parameters ---------------------------------------------
tFinal  = 1e-2;
dtList  = [5e-6 2e-7];          % coarse → finest (all ultra-small)
boundary = "neumann";

% SAV parameters
gamma0 = 4;  eta = 0.95;  xi_flag = 1;  C0 = 0;

Ndisp = numel(dtList);               % show every step

%% 3.  Storage ------------------------------------------------------------
n = numel(dtList);
phiNMG = cell(1,n);  phiSAV = cell(1,n);

%% 4.  Run solvers --------------------------------------------------------
for k = 1:n
    dt  = dtList(k);
    nSt = ceil(tFinal/dt);
    
    % --- NMG
    [~,tmp,~,~] = CahnHilliard_NMG(phi0,'t_iter',nSt,'dt',dt,...
                                   'dt_out',nSt,'boundary',boundary);
    phiNMG{k} = tmp(:,:,end);
    
    % --- SAV
    [~,tmp,~,~] = CahnHilliard_SAV(phi0,'t_iter',nSt,'dt',dt,...
                                   'dt_out',nSt,'boundary',boundary,...
                                   'gamma0',gamma0,'eta',eta,...
                                   'xi_flag',xi_flag,'C0',C0);
    phiSAV{k} = tmp(:,:,end);
    
    fprintf('done Δt = %.1e\n',dt);
end

%% 5.  Build triplet storyboard ------------------------------------------
rbMap = interp1(1:100:1100, redbluecmap, 1:1001);   % diverging map
cLim  = [-1 1];                                     % phase-field range

% common limits for differences
maxDiff = 0;
for k = 1:n
    maxDiff = max(maxDiff, max(abs(phiNMG{k}(:)-phiSAV{k}(:))));
end
cLimD = maxDiff * [-1 1];

figure('Name','NMG | SAV | difference (ultra-fine steps)', ...
       'Position',[80 80 900 280*Ndisp]);

tiledlayout(Ndisp,3,'Padding','compact','TileSpacing','compact');

for r = 1:Ndisp
    dt = dtList(r);
    nF = phiNMG{r};  sF = phiSAV{r};  dF = nF - sF;
    
    % NMG field
    nexttile;
    imagesc(nF); axis equal tight off;
    colormap(gca,rbMap); caxis(cLim);
    if r==1, title('NMG'); end
    ylabel(sprintf('\\Delta t = %.1e',dt));
    
    % SAV field
    nexttile;
    imagesc(sF); axis equal tight off;
    colormap(gca,rbMap); caxis(cLim);
    if r==1, title('SAV'); end
    
    % Difference
    nexttile;
    imagesc(dF); axis equal tight off;
    colormap(gca,'jet'); caxis(cLimD);
    if r==1, title('\phi_{NMG} - \phi_{SAV}'); end
end

% colour bars
cb1 = colorbar('Position',[0.92 0.57 0.015 0.35]); cb1.Label.String='\phi';
cb2 = colorbar('Position',[0.92 0.12 0.015 0.35]); cb2.Label.String='\Delta\phi';

% save optional
% saveas(gcf,'CH_triplet_ultrafine.png');
% FIGURE 1
indir = "../IC/";
outdir = "../output/output_MATLAB-periodic";

m = 8;
GridSize = 64;
h = 1/GridSize;
epsilon = m * h/ (2 * sqrt(2) * atanh(0.9));
dt = 5.5e-6;
max_it = 10;
boundary = 'periodic';

print_phi = true;
dt_out = 1;
nx = GridSize;
ny = GridSize;

phi0 = zeros(nx,ny);

phi0(1,10) = 2;
phi0(64,25) = 4;
phi0(32, 1) = 6;
phi0(50, 64) = 8;
phi0(32,32) = 10;
                
% % #################################################
% % RUN NMG SOLVER 
% % #################################################

% pathname = sprintf("%s/periodic_NMG_test_",outdir);
% fprintf("Running NMG solver with parameters: %s\n", pathname);
% tStart_NMG = tic;
% [t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_NMG(phi0,...
%                                     t_iter = max_it,...
%                                     dt = dt,...
%                                     m = m,...
%                                     boundary = boundary,...
%                                     printphi=print_phi,...
%                                     pathname=pathname,...
%                                     dt_out = dt_out);
% elapsedTime = toc(tStart_NMG);


% writematrix(delta_mass_t,sprintf('%smass_uncentered.csv', pathname));
% writematrix(E_t,sprintf('%senergy.csv', pathname));

pathname = sprintf("%s/periodic_SAV_test_",outdir);
fprintf("Running SAV solver with parameters: %s\n", pathname);
tStart_SAV = tic;
[t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_SAV(phi0,...
                                    t_iter = max_it,...
                                    dt = dt,...
                                    m = m,...
                                    boundary = boundary,...
                                    printphi=print_phi,...
                                    pathname=pathname,...
                                    dt_out = dt_out);
elapsedTime = toc(tStart_SAV);

writematrix(delta_mass_t,sprintf('%smass_uncentered.csv', pathname));
writematrix(E_t,sprintf('%senergy.csv', pathname));

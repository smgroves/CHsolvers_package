% FIGURE 1
indir = "../IC/";
boundary = 'periodic';

outdir = sprintf("../output/output_MATLAB-%s", boundary);

m = 8;
GridSize = 128;
h = 1/GridSize;
epsilon = m * h/ (2 * sqrt(2) * atanh(0.9));
dt = 5.5e-6;
max_it = 200;

print_phi = true;
dt_out = 1;
nx = GridSize;
ny = GridSize;
n_relax = 4;
% #################################################
init_file = sprintf("%s/initial_phi_%d_smooth_n_relax_%d.csv",indir,GridSize, n_relax);
phi0 = readmatrix(init_file);

% % #################################################
% % RUN NMG SOLVER 
% % #################################################

pathname = sprintf("%s/%s_NMG_finaltest_",outdir, boundary);
fprintf("Running NMG solver with parameters: %s\n", pathname);
tStart_NMG = tic;
[t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_NMG(phi0,...
                                    t_iter = max_it,...
                                    dt = dt,...
                                    m = m,...
                                    boundary = boundary,...
                                    printphi=print_phi,...
                                    pathname=pathname,...
                                    dt_out = dt_out);
elapsedTime = toc(tStart_NMG);


writematrix(delta_mass_t,sprintf('%smass.csv', pathname));
writematrix(E_t,sprintf('%senergy.csv', pathname));


% % #################################################
% % RUN SAV SOLVER 
% % #################################################

pathname = sprintf("%s/%s_SAV_finaltest_",outdir, boundary);
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

writematrix(delta_mass_t,sprintf('%smass.csv', pathname));
writematrix(E_t,sprintf('%senergy.csv', pathname));

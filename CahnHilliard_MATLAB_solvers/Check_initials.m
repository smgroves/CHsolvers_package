%% Define file paths for the initial conditions
file_path128 = '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_128_smooth_n_relax_4.csv';
file_path256 = '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_256_smooth_n_relax_4.csv';
file_path512 = '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_512_smooth_n_relax_4.csv';

%% Load the initial condition files (assumed to be 2D arrays)
phi128 = readmatrix(file_path128);  % size: 128 x 128
phi256 = readmatrix(file_path256);  % size: 256 x 256
phi512 = readmatrix(file_path512);  % size: 512 x 512

%% Determine the common grid points
% For periodic BCs over [0,1), the grid points are defined as:
%   x = 0, 1/N, 2/N, ..., (N-1)/N.
% Thus, the 128-grid is a subset of the 256-grid (every 2nd point)
% and the 512-grid (every 4th point).

% Extract the common subset from phi256 and phi512:
phi256_common = phi256(1:2:end, 1:2:end);  % 256 -> 128
phi512_common = phi512(1:4:end, 1:4:end);  % 512 -> 128

%% Compute the relative L2 differences at the common grid points
% Using the 128-grid as the reference:
rel_diff_128_256 = norm(phi128(:) - phi256_common(:), 2) / norm(phi128(:),2);
rel_diff_128_512 = norm(phi128(:) - phi512_common(:), 2) / norm(phi128(:),2);

fprintf('Relative L2 difference between 128 and 256 initial conditions: %.3e\n', rel_diff_128_256);
fprintf('Relative L2 difference between 128 and 512 initial conditions: %.3e\n', rel_diff_128_512);

%% Visualize the common points for visual inspection
% Create common grids for plotting (for a periodic domain [0,1))
N_common = 128;
x_common = linspace(0, 1, N_common+1); x_common = x_common(1:end-1);
y_common = linspace(0, 1, N_common+1); y_common = y_common(1:end-1);
[X_common, Y_common] = meshgrid(x_common, y_common);

figure;
subplot(1,3,1);
imagesc(x_common, y_common, phi128);
axis equal tight; colorbar;
title('Initial Condition: 128x128');

subplot(1,3,2);
imagesc(x_common, y_common, phi256_common);
axis equal tight; colorbar;
title('256x256 (common points)');

subplot(1,3,3);
imagesc(x_common, y_common, phi512_common);
axis equal tight; colorbar;
title('512x512 (common points)');

sgtitle('Comparison of 2D Initial Conditions at Common Grid Points');
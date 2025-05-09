%% Plotting Unit-Circle Initial Conditions on Nested Grids
clear; close all; clc;

% Domain and initial condition parameters
Lx = 1; Ly = 1;
R0 = 0.1;      % Interface radius
m = 8;         % Parameter (as before)

% Resolutions to test
N_list = [128, 256, 512];
num_res = length(N_list);

% Preallocate cell arrays for the initial conditions and grids
phi0_all = cell(num_res, 1);
grids = cell(num_res, 1);

for i = 1:num_res
    % Set grid resolution
    Nx = N_list(i);
    Ny = Nx;
    hx = Lx / Nx;
    hy = Ly / Ny;
    
    % Create periodic grid: points in [0,1)
    x = linspace(0, 1, Nx+1); x = x(1:end-1);
    y = linspace(0, 1, Ny+1); y = y(1:end-1);
    [X, Y] = meshgrid(x, y);
    grids{i} = struct('X', X, 'Y', Y);
    
    % Compute gam and interface thickness delta
    gam = m * hx / (2 * sqrt(2) * atanh(0.9));
    eps_c = gam;
    delta = eps_c * sqrt(2);
    
    % Compute distance from the center (0.5, 0.5)
    R = sqrt((X - 0.5).^2 + (Y - 0.5).^2);
    
    % Compute the tanh profile for psi0 and then phi0
    psi0 = 0.5 * (1 - tanh((R - R0) / (2 * delta)));
    phi0 = 2 * psi0 - 1;
    
    % Store the computed initial condition
    phi0_all{i} = phi0;
    
    fprintf('Computed initial condition for Nx = %d, hx = %.3e, gam = %.3e\n', Nx, hx, gam);
end

%% Plot the three initial conditions side by side for visual checking
figure;
for i = 1:num_res
    subplot(1, num_res, i);
    imagesc(grids{i}.X(1,:), grids{i}.Y(:,1), phi0_all{i});
    axis equal tight; colorbar;
    title(sprintf('Nx = %d', N_list(i)));
end
sgtitle('Unit Circle Initial Conditions for Different Resolutions');


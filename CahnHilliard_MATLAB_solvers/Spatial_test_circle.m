%% Spatial Convergence Test Using SAV Solver with a Fixed Continuous Initial Condition
clear; close all; clc;

%% Domain and Continuous Initial Condition Parameters
Lx = 1; Ly = 1;       % Domain size
R0 = 0.1;             % Fixed interface radius
delta = 0.01;         % Fixed interface thickness

% Define the continuous initial condition function (unit circle with tanh profile)
% psi0 = 0.5*(1 - tanh((R - R0)/(2*delta)))  and  phi0 = 2*psi0 - 1.
phi0_fun = @(X,Y) 2*(0.5*(1 - tanh((sqrt((X-0.5).^2 + (Y-0.5).^2) - R0) / (2*delta)))) - 1;

%% SAV Solver Parameters
total_time = 2e-1;    % Total simulation time
dt = 1e-4;            % Time step
dt_out = 10;          % Output frequency (choose a reasonable value to record outputs)
m_SAV = 8;            % Parameter for the SAV solver
Beta = 0;             % Relaxation parameter (as before)
boundaryChoice = 'periodic';

%% Resolutions for the Spatial Convergence Test
N_list = [128, 256, 512];  % Grid resolutions (Nx = Ny)
num_res = length(N_list);

% Preallocate arrays/cells for storing final solutions and grid info
phi_final = cell(num_res, 1);  % Final solution from SAV solver
grids = cell(num_res, 1);      % Grid structure (X and Y)
Nx_values = zeros(num_res, 1); % Number of grid points

%% Loop Over Each Resolution: Build Initial Condition, Run SAV Solver, Store Final Solution
for i = 1:num_res
    % Set grid resolution
    Nx = N_list(i);
    Ny = Nx;
    Nx_values(i) = Nx;
    hx = Lx / Nx;
    hy = Ly / Ny;
    
    % Create grid for a periodic domain [0,1) in both x and y.
    x = linspace(0, 1, Nx+1); x = x(1:end-1);
    y = linspace(0, 1, Ny+1); y = y(1:end-1);
    [X, Y] = meshgrid(x, y);
    grids{i} = struct('X', X, 'Y', Y);
    
    % Sample the continuous initial condition on the grid
    phi0 = phi0_fun(X, Y);
    
    % Run the SAV solver using this initial condition.
    max_it = round(total_time / dt);
    % Call the SAV solver (assumed to be in your path)
    [t_out, phi_t, ~, ~, ~] = CahnHilliard_SAV(phi0, ...
                                    't_iter', max_it, ...
                                    'dt', dt, ...
                                    'm', m_SAV, ...
                                    'dt_out', dt_out, ...
                                    'boundary', boundaryChoice, ...
                                    'Beta', Beta);
    % Assume phi_t is a 3D array (Nx x Ny x time); take the final time slice.
    phi_final{i} = phi_t(:, :, end);
    
    fprintf('Resolution: %dx%d, hx = %.3e\n', Nx, Ny, hx);
end

%% Use the Finest Resolution as the Reference Solution
[~, idx_ref] = max(Nx_values);
phi_ref = phi_final{idx_ref};
grid_ref = grids{idx_ref};
Nx_ref = Nx_values(idx_ref);

%% Compute Relative L2 Errors at Common Grid Points (Nested Grids)
% Since the grids are nested (e.g. 128x128 is contained in 512x512),
% we extract the corresponding points from the finest-grid solution.
errors = zeros(num_res, 1);
for i = 1:num_res
    if i == idx_ref
        errors(i) = NaN;  % Skip the reference solution
        continue;
    end
    
    % Extraction factor: how many fine-grid points per coarse-grid interval.
    factor = Nx_ref / Nx_values(i);
    if mod(factor, 1) ~= 0
        error('Finest grid size is not an integer multiple of the current grid size.');
    end
    
    % Extract common grid points from the reference solution.
    phi_ref_common = phi_ref(1:factor:end, 1:factor:end);
    
    % Compute the relative L2 error norm:
    % relative error = norm(coarse - common_ref) / norm(coarse)
    rel_err = norm(phi_final{i}(:) - phi_ref_common(:), 2) / norm(phi_final{i}(:), 2);
    errors(i) = rel_err;
    
    fprintf('Nx = %d, extraction factor = %d, relative L2 error = %.3e\n', ...
        Nx_values(i), factor, rel_err);
end

%% Plot Spatial Convergence: Error vs. Grid Spacing
% For a periodic grid, grid spacing h = 1/Nx.
hs = 1 ./ Nx_values;
valid_idx = ~isnan(errors);

figure;
loglog(hs(valid_idx), errors(valid_idx), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Grid spacing, h');
ylabel('Relative L_2 Error Norm');
title(sprintf('Spatial Convergence Test using SAV Solver (dt = %.1e)', dt));
grid on;
hold on;

% Linear fit in log-log space to estimate the convergence order
p = polyfit(log(hs(valid_idx)), log(errors(valid_idx)), 1);
fitted_line = exp(polyval(p, log(hs(valid_idx))));
loglog(hs(valid_idx), fitted_line, '--', 'LineWidth', 2);
legend('Measured Error', sprintf('Fit: Slope = %.2f', p(1)), 'Location', 'northeastoutside');

fprintf('Estimated spatial convergence order (slope): %.2f\n', p(1));

%% Optionally: Plot Initial and Final Solutions for Visual Inspection
figure;
for i = 1:num_res
    % Plot the initial condition on the top row
    subplot(2, num_res, i);
    imagesc(grids{i}.X(1,:), grids{i}.Y(:,1), phi0_fun(grids{i}.X, grids{i}.Y));
    axis equal tight; colorbar;
    title(sprintf('Initial, Nx=%d', N_list(i)));
    
    % Plot the final SAV solution on the bottom row
    subplot(2, num_res, i + num_res);
    imagesc(grids{i}.X(1,:), grids{i}.Y(:,1), phi_final{i});
    axis equal tight; colorbar;
    title(sprintf('Final, Nx=%d', N_list(i)));
end
sgtitle('Initial and Final Solutions from SAV Solver');
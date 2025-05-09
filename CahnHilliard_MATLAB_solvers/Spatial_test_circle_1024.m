%% Spatial Convergence Test Using SAV Solver with Extended Resolutions
clear; close all; clc;

%% Domain and Continuous Initial Condition Parameters
Lx = 1; Ly = 1;       % Domain size
R0 = 0.1;             % Fixed interface radius
delta = 0.01;         % Fixed interface thickness

% Define the continuous initial condition (unit circle with tanh profile)
phi0_fun = @(X,Y) 2*(0.5*(1 - tanh((sqrt((X-0.5).^2 + (Y-0.5).^2) - R0) / (2*delta)))) - 1;

%% SAV Solver Parameters
total_time = 2e-2;    % Total simulation time
dt = 1e-4;            % Time step
dt_out = 10;          % Output frequency (adjust as needed)
m_SAV = 8;            % Parameter for the SAV solver
Beta = 0;             % Relaxation parameter
boundaryChoice = 'periodic';

%% Resolutions for the Spatial Convergence Test
N_list = [128, 256, 512, 1024];  % Extended list of grid resolutions
num_res = length(N_list);

% Preallocate cells for final solutions and grid data
phi_final = cell(num_res, 1);   % Final SAV solutions
grids     = cell(num_res, 1);   % Grid coordinates for each resolution
Nx_values = zeros(num_res, 1);   % Number of grid points (assumed square)

%% Loop Over Each Resolution: Generate IC, Run SAV Solver, Store Final Solution
for i = 1:num_res
    % Set grid resolution
    Nx = N_list(i);
    Ny = Nx;
    Nx_values(i) = Nx;
    hx = Lx / Nx;
    hy = Ly / Ny;
    
    % Create grid for periodic domain [0,1) in both x and y
    x = linspace(0, 1, Nx+1); x = x(1:end-1);
    y = linspace(0, 1, Ny+1); y = y(1:end-1);
    [X, Y] = meshgrid(x, y);
    grids{i} = struct('X', X, 'Y', Y);
    
    % Sample the fixed continuous initial condition on the grid
    phi0 = phi0_fun(X, Y);
    
    % Run the SAV solver with this initial condition
    max_it = round(total_time / dt);
    [t_out, phi_t, ~, ~, ~] = CahnHilliard_SAV(phi0, ...
                                    't_iter', max_it, ...
                                    'dt', dt, ...
                                    'm', m_SAV, ...
                                    'dt_out', dt_out, ...
                                    'boundary', boundaryChoice, ...
                                    'Beta', Beta);
    % Extract the final solution (last time slice)
    phi_final{i} = phi_t(:, :, end);
    
    fprintf('Resolution: %dx%d, hx = %.3e\n', Nx, Ny, hx);
end

%% Select the Finest Resolution (1024x1024) as the Reference Solution
[~, idx_ref] = max(Nx_values);
phi_ref = phi_final{idx_ref};
grid_ref = grids{idx_ref};
Nx_ref = Nx_values(idx_ref);

%% Compute Relative L2 Errors Using Nested Extraction
errors = zeros(num_res, 1);
for i = 1:num_res
    if i == idx_ref
        errors(i) = NaN;  % Skip comparison for the reference solution
        continue;
    end
    
    % Compute extraction factor: how many fine-grid points per coarse-grid interval
    factor = Nx_ref / Nx_values(i);
    if mod(factor, 1) ~= 0
        error('Finest grid size is not an integer multiple of the current grid size.');
    end
    
    % Extract common grid points from the reference solution
    phi_ref_common = phi_ref(1:factor:end, 1:factor:end);
    
    % Compute the relative L2 error norm: norm(coarse - common_ref)/norm(coarse)
    rel_err = norm(phi_final{i}(:) - phi_ref_common(:), 2) / norm(phi_final{i}(:), 2);
    errors(i) = rel_err;
    
    fprintf('Nx = %d, extraction factor = %d, relative L2 error = %.3e\n', ...
        Nx_values(i), factor, rel_err);
end

%% Plot Spatial Convergence: Error vs. Grid Spacing
hs = 1 ./ Nx_values;  % For periodic grid, h = 1/N
valid_idx = ~isnan(errors);

figure;
loglog(hs(valid_idx), errors(valid_idx), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Grid spacing, h');
ylabel('Relative L_2 Error Norm');
title(sprintf('Spatial Convergence Test (dt = %.1e)', dt));
grid on;
hold on;

% Linear fit in log-log space to estimate convergence order
p = polyfit(log(hs(valid_idx)), log(errors(valid_idx)), 1);
fitted_line = exp(polyval(p, log(hs(valid_idx))));
loglog(hs(valid_idx), fitted_line, '--', 'LineWidth', 2);
legend('Measured Error', sprintf('Fit: Slope = %.2f', p(1)), 'Location', 'northeastoutside');

fprintf('Estimated spatial convergence order (slope): %.2f\n', p(1));

%% Optional: Plot Initial and Final Solutions for Visual Inspection
figure;
for i = 1:num_res
    subplot(2, num_res, i);
    imagesc(grids{i}.X(1,:), grids{i}.Y(:,1), phi0_fun(grids{i}.X, grids{i}.Y));
    axis equal tight; colorbar;
    title(sprintf('IC, Nx=%d', N_list(i)));
    
    subplot(2, num_res, i+num_res);
    imagesc(grids{i}.X(1,:), grids{i}.Y(:,1), phi_final{i});
    axis equal tight; colorbar;
    title(sprintf('Final, Nx=%d', N_list(i)));
end
sgtitle('Initial and Final Solutions from SAV Solver');
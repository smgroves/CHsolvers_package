%% Spatial Convergence Test for 2D Data (Periodic BCs) Using Nested Grids
% Define file paths for the initial condition files (assumed to be 128x128, 256x256, 512x512)
file_paths = { ...
    '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_128_smooth_n_relax_4.csv', ...
    '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_256_smooth_n_relax_4.csv', ...
    '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_512_smooth_n_relax_4.csv'};

num_files = length(file_paths);
phi_final = cell(num_files,1);
Nx_values = zeros(num_files,1);
grids = cell(num_files,1);

%% Simulation Parameters
total_time = 2e-1;   % Total simulation time
dt         = 1e-4;   % Time step (same for all resolutions)
dt_out     = 1e+3;   % Output frequency for the solver
m          = 8;      % Parameter for the SAV solver
boundaryChoice = 'periodic';
Beta = 0e3;

%% Loop over each resolution: load initial data, run simulation, and store final solution
for i = 1:num_files
    % Load the 2D initial condition
    phi0 = readmatrix(file_paths{i});
    [Nx, Ny] = size(phi0);
    Nx_values(i) = Nx;
    
    % Define grid for periodic domain [0,1)^2:
    % Use N points spanning [0,1) in each direction
    x = linspace(0,1,Nx+1); x = x(1:end-1);
    y = linspace(0,1,Ny+1); y = y(1:end-1);
    [X, Y] = meshgrid(x,y);
    grids{i} = struct('X', X, 'Y', Y);
    
    % Run the SAV simulation; assume the final solution is the last time slice of phi_t
    max_it = round(total_time/dt);
    [~, phi_t, ~, ~, ~] = CahnHilliard_SAV(phi0, ...
                              't_iter', max_it, ...
                              'dt', dt, ...
                              'm', m, ...
                              'dt_out', dt_out, ...
                              'boundary', boundaryChoice, ...
                              'Beta', Beta);
    phi_final{i} = phi_t(:,:,end);
    
    fprintf('Loaded file: %s with Nx = %d, Ny = %d\n', file_paths{i}, Nx, Ny);
end

%% Choose the finest resolution (largest Nx) as the reference solution
[~, idx_ref] = max(Nx_values);
phi_ref = phi_final{idx_ref};
grid_ref = grids{idx_ref};
Nx_ref = Nx_values(idx_ref);

%% Compute errors for spatial convergence using nested extraction
% For periodic BCs, grid spacing is defined as h = 1/Nx.
errors = zeros(num_files, 1);
hs = zeros(num_files, 1);

for i = 1:num_files
    Nx_i = Nx_values(i);
    h = 1 / Nx_i;  % grid spacing for the current resolution
    hs(i) = h;
    
    if i == idx_ref
        errors(i) = NaN;  % Skip error computation for the reference case
        continue;
    end
    
    % Determine the extraction factor.
    % For example, if the reference grid is 512x512 and current grid is 256x256,
    % then factor = 512/256 = 2; for 128x128, factor = 512/128 = 4.
    factor = Nx_ref / Nx_i;
    if mod(factor, 1) ~= 0
        error('Finest grid size is not an integer multiple of the current grid size.');
    end
    
    % Extract the common grid points from the reference solution.
    % This gives a "reference" solution on the same grid as the current resolution.
    phi_ref_common = phi_ref(1:factor:end, 1:factor:end);
    
    % Compute the L2 error norm over the domain.
    % Multiply by h (area element) to approximate the continuous L2 norm.
    error_i = sqrt(sum(sum((phi_final{i} - phi_ref_common).^2)) * h^2);
    errors(i) = error_i;
    
    fprintf('Nx = %d, h = %.3e, L2 Error = %.3e\n', Nx_i, h, error_i);
end

%% Remove the reference (NaN) from errors and hs for plotting
valid_idx = ~isnan(errors);
hs_plot = hs(valid_idx);
errors_plot = errors(valid_idx);

%% Plot spatial convergence: error vs grid spacing on a log-log plot
figure;
loglog(hs_plot, errors_plot, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Grid spacing, h');
ylabel('L_2 Error Norm');
title(sprintf('Spatial Convergence Test (%s)\nTemporal dt = %.1e', boundaryChoice, dt));
grid on;
hold on;

% Perform a linear fit in log-log space to estimate the spatial order
p = polyfit(log(hs_plot), log(errors_plot), 1);
fitted_line = exp(polyval(p, log(hs_plot)));

% Plot the fitted line and place the legend outside the plot area
loglog(hs_plot, fitted_line, '--', 'LineWidth', 2);
legend('Measured Error', ['Fit: Slope = ' num2str(p(1), '%.2f')], 'Location', 'northeastoutside');

fprintf('Estimated Spatial Order of Convergence (slope): %.2f\n', p(1));
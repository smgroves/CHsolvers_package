%% Stability Test: dt Sensitivity via D (Auxiliary Variable Error)
clear; close all; clc;

%% Domain & Initial Condition 

% Define the continuous initial condition (circle profile)
% R0 = 0.1; 
% delta = 0.01;  % Fixed interface thickness
% phi0_fun = @(X,Y) 2*(0.5*(1 - tanh((sqrt((X-0.5).^2+(Y-0.5).^2)-R0)/(2*delta)))) - 1;
% 
% N = 128;
% x = linspace(0,1,N+1); x = x(1:end-1);  % Periodic grid in x
% y = linspace(0,1,N+1); y = y(1:end-1);  % Periodic grid in y
% [X,Y] = meshgrid(x,y);
% phi0 = phi0_fun(X,Y);

% Define the continuous initial condition (spinodal profile)
file_path = '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_128_smooth_n_relax_4.csv';
phi0 = readmatrix(file_path);

%% SAV Solver Parameters for Stability Test
total_time = 2e-2;         % Total simulation time
dt_ref = 5e-6;             % Reference smallest time step
dt_list = dt_ref * [2, 4, 8];  % dt_list: [1e-5, 2e-5, 4e-5]
dt_out  = 1e2;             % Output every time step
m = 8;                       % SAV parameter
Beta = 0e3;                  % Relaxation parameter
eta = 0;                     % Relaxation parameter
boundaryChoice = 'neumann';      

num_dt = length(dt_list);
max_D = zeros(num_dt, 1);
D_all = cell(num_dt, 1);
t_all = cell(num_dt, 1);
phi_final_all = cell(num_dt, 1);  % To store final phi for each dt

%% Loop over different dt values and run the SAV solver
for k = 1:num_dt
    dt = dt_list(k);
    t_iter = round(total_time/dt);
    
    [t_out, phi_t, ~, ~, D_t] = CahnHilliard_SAV(phi0, ...
         't_iter', t_iter, ...
         'dt', dt, ...
         'dt_out', dt_out, ...
         'm', m, ...
         'boundary', boundaryChoice, ...
         'eta', eta);
     
    max_D(k) = max(abs(D_t));
    D_all{k} = D_t;
    t_all{k} = t_out;
    phi_final_all{k} = phi_t(:,:,end);  % Store final phi
    
    fprintf('dt = %.1e: max D = %.3e\n', dt, max(abs(D_t)));
end

%% Plot Stability: Maximum D vs. dt
figure(1);
loglog(dt_list, max_D, 's-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', sprintf('eta = %.1e', eta));
xlabel('Time step, dt');
ylabel('Maximum |D|');
title('Relaxed_SAV Test: Maximum D vs. dt');
grid on;

%% Plot Time Evolution of D for All dt
figure(2);
hold on;
colors = lines(num_dt);
for k = 1:num_dt
    plot(t_all{k}, D_all{k}, 'Color', colors(k,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('eta = %.1e, dt = %.1e', eta, dt_list(k)));
end
xlabel('Time');
ylabel('D (Auxiliary Variable Error)');
title('Time Evolution of D for Different dt');
legend('Location', 'best');
grid on;
hold off;

%% Plot Final φ for All dt
% We assume that the final phi is on the same physical grid.
% Create a grid from the CSV (N=128)
[Nx, Ny] = size(phi0);
x = linspace(0, 1, Nx+1); x = x(1:end-1);
y = linspace(0, 1, Ny+1); y = y(1:end-1);

figure(3);
for k = 1:num_dt
    subplot(1, num_dt, k);
    imagesc(x, y, phi_final_all{k});
    axis equal tight; colorbar;
    title(sprintf('Final \\phi, dt = %.1e', dt_list(k)));
end
sgtitle('Final \phi for Different dt (CSV Initial, N=128)');
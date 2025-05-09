%% Set up the initial parameters
% File with initial condition
    file_path = '/Users/luminjhe/Downloads/CHsolvers_package/IC/initial_phi_128_smooth_n_relax_4.csv';
    phi0 = readmatrix(file_path);
    
    % total_time = 1e-1;   % Total simulation time
    % dt_out     = 1e+3;   % Output frequency (as in your solver)
    total_time = 2e-2;   % Total simulation time
    dt_out     = 1e+1;   % Output frequency (as in your solver)
    boundaryChoice = 'neumann';  % Pre-determine the boundary condition

%% Compute the reference solution using the smallest dt
    dt_ref = 1e-6;
    max_it_ref = round(total_time / dt_ref);
    [t_out_ref, phi_t_ref, delta_mass_t_ref, E_t_ref, D_t_ref] = ...
        CahnHilliard_SAV(phi0, ...
                         't_iter', max_it_ref, ...
                         'dt', dt_ref, ...
                         'm', 8, ...
                         'dt_out', dt_out, ...
                         'boundary', boundaryChoice, ...
                         'Beta', 1e3);
                 
% Assume the solver returns phi_t_ref as a matrix where each column corresponds
% to the solution at a given output time. Take here the final time column.
    phi_ref_final = phi_t_ref(:, end);

%% Define a set of coarser time steps using a geometric progression
    dt_list = dt_ref * [2, 4, 8];
    num_tests = length(dt_list);
    errors = zeros(num_tests, 1);

%% Loop over each dt to compute the error relative to the reference solution
    for i = 1:num_tests
        dt_current = dt_list(i);
        max_it_current = round(total_time / dt_current);
        
        [~, phi_t_current, ~, ~, ~] = ...
            CahnHilliard_SAV(phi0, ...
                             't_iter', max_it_current, ...
                             'dt', dt_current, ...
                             'm', 8, ...
                             'dt_out', dt_out, ...
                             'boundary', boundaryChoice, ...
                             'Beta', 0e3);
        
        % Get the final time solution for the current dt
        phi_current_final = phi_t_current(:, end);
        
        % Compute the L2 norm of the difference between the current and reference solutions
        errors(i) = norm(phi_current_final - phi_ref_final, 2);
        
        fprintf('dt = %.1e, L2 error = %.3e\n', dt_current, errors(i));
    end

%% Plot the convergence results on a log-log scale
    
    % Set additional parameters for the title
    Nx = size(phi0, 1);  % Number of spatial grid points (assuming phi0 is a column vector)
    finest_dt = dt_ref;       % Finest time step used
    
    figure;
    loglog(dt_list, errors, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Time step, dt');
    ylabel('L_2 Error Norm');
    
    % Create a two-line title using \n to break the line
    title(sprintf('Temporal Convergence\n(%s, Nx = %d)\n(Finest dt = %.1e)', ...
          boundaryChoice, Nx, finest_dt));
    
    grid on;
    hold on;
    
    % Perform a linear fit in log-log space
    p = polyfit(log(dt_list), log(errors), 1);
    fitted_line = exp(polyval(p, log(dt_list)));
    
    % Plot the fitted line
    loglog(dt_list, fitted_line, '--', 'LineWidth', 2);
    
    % Annotate the slope (which is the experimental order of convergence)
    legend('Measured Error', ['Fit: Slope = ' num2str(p(1), '%.2f')], ...
           'Location', 'northeastoutside');
    
    % Display the estimated order of convergence in the console
    fprintf('Estimated Order of Convergence (slope): %.2f\n', p(1));

%% Compute the Experimental Order of Convergence (EOC)
    EOC = log(errors(1:end-1)./errors(2:end)) ./ log(dt_list(1:end-1)./dt_list(2:end));
    disp('Experimental Order of Convergence (EOC) for each dt reduction:');
    disp(EOC);


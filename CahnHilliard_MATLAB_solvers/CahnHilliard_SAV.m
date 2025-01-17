function [t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_SAV(phi0, varargin)
% This function uses the scalar auxiliary variable method to solve the 
% Cahn-Hilliard equation for a specified number of time steps of size dt.
% 
% INPUTS
    % phi0 = Initial field of chemical states in the domain, created by ch_initialization.
%
%NAME-VALUE PAIRS
    % t_iter = Number of time steps simulated. Default = 1e3.
    % dt = Time step. Default = 2.5e-5 characteristic times.
    % dt_out = Every 'dt_out' time step data is stored. Default = 10.
    % m = Number of mesh points over which the interface exists. Default = 4.
    % epsilon2 = Squared interfacial transition distance; if specified,
    %   m will be overwritten. Default = nan (do not overwrite m).
    % boundary = Boundary conditions for the simulation:
    %   'periodic' (default) - flux on one domain border equals negative flux on the opposite border.
    %   'neumann' - zero flux on the domain borders.
    % domain = Vector of rightmost and leftmost grid points in x and y.
    %   Format: [xright xleft yright yleft]. Default = [1 0 1 0].
    % printphi = Logical to print phi to a file. Default = true.
    % pathname = Name of the path to which phi is printed. Default = 'cd'.
%
%OUTPUT
    % t_out = Time corresponding to the dt time step outputs.
    % phi_t = Multidimensional array of phi over t_out.
    % delta_mass_t = Vector of mass change over t_out.
    % E_t = Vector of total energy over t_out.
    
%% Set option defaults and parse inputs

    % Set parameter defaults
    default_t_iter = 1e3;
    default_dt = 2.5e-5;
    default_dt_out = 10;
    default_m = 4;
    default_epsilon2 = nan;
    default_boundary = 'periodic';
    default_domain = [1 0 1 0];
    default_printphi = true;
    default_pathname = 'cd';
    default_C0 = 0;
    default_Beta = 0;
    default_gamma0 = 0;
    
    CahnHilliard_SAV_parser = inputParser;

    % Set general criteria for inputs and name-value pairs
    valid_matrix = @(x) ismatrix(x);
    valid_integer = @(x) x-floor(x) == 0;
    valid_pos_num = @(x) isnumeric(x) && (x > 0);
    valid_boundary_type = @(x) strcmpi(x,'periodic') || strcmpi(x,'neumann');
    valid_domain_vector = @(x) length(x) == 4;
    valid_logical = @(x) islogical(x) || x == 1 || x == 0;
    valid_string = @(x) ischar(x) || isstring(x);


    % Set parser options and valid input criteria
    addRequired(CahnHilliard_SAV_parser,'phi0',valid_matrix);
    
    addParameter(CahnHilliard_SAV_parser,'t_iter',default_t_iter,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'dt',default_dt,valid_pos_num);
    addParameter(CahnHilliard_SAV_parser,'dt_out',default_dt_out,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'m',default_m,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'epsilon2',default_epsilon2,valid_pos_num);
    addParameter(CahnHilliard_SAV_parser,'domain',default_domain,valid_domain_vector);
    addParameter(CahnHilliard_SAV_parser,'boundary',default_boundary,valid_boundary_type);
    addParameter(CahnHilliard_SAV_parser,'printphi',default_printphi,valid_logical);
    addParameter(CahnHilliard_SAV_parser,'pathname',default_pathname,valid_string);
    addParameter(CahnHilliard_SAV_parser,'C0',default_C0,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'Beta',default_Beta,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'gamma0',default_gamma0,valid_integer);

    parse(CahnHilliard_SAV_parser, phi0, varargin{:});
    
    % Extract parsed inputs
    phi0 = CahnHilliard_SAV_parser.Results.phi0;
    t_iter = CahnHilliard_SAV_parser.Results.t_iter;
    dt = CahnHilliard_SAV_parser.Results.dt;
    dt_out = CahnHilliard_SAV_parser.Results.dt_out;
    if dt_out > t_iter, error('Error: dt_out should not be greater than t_iter.'); end
    m = CahnHilliard_SAV_parser.Results.m;
    epsilon2 = CahnHilliard_SAV_parser.Results.epsilon2;
    boundary = CahnHilliard_SAV_parser.Results.boundary;
    xright = CahnHilliard_SAV_parser.Results.domain(1);
    xleft = CahnHilliard_SAV_parser.Results.domain(2);
    yright = CahnHilliard_SAV_parser.Results.domain(3);
    yleft = CahnHilliard_SAV_parser.Results.domain(4);
    printphi = CahnHilliard_SAV_parser.Results.printphi;
    pathname = CahnHilliard_SAV_parser.Results.pathname;
    C0 = CahnHilliard_SAV_parser.Results.C0;
    Beta = CahnHilliard_SAV_parser.Results.Beta;
    gamma0 = CahnHilliard_SAV_parser.Results.gamma0;

%% Define and initialize key simulation parameters

    [nx,ny] = size(phi0); % Define number of grid points in x and y
    Lx = xright-xleft; Ly = yright-yleft;
    hx = Lx/nx; hy = Ly/ny;
    h2 = hx*hy; % Define mesh size
    k_x = 1i*[0:nx/2 -nx/2+1:-1]*(2*pi/Lx); k_y = 1i*[0:ny/2 -ny/2+1:-1]*(2*pi/Ly);
    k_xx = k_x.^2; k_yy = k_y.^2;
    [kxx,kyy] = meshgrid(k_xx,k_yy);
    k2 = kxx + kyy;
    k4 = k2.^2;
    if isnan(epsilon2)
        epsilon2 = h2*m^2/(2*sqrt(2)*atanh(0.9))^2; % Define ϵ^2 if not prespecified
    else
        m = sqrt((epsilon2*(2*sqrt(2)*atanh(0.9))^2)/h2); % Else overwrite m
        display(m);
    end

%% Initialization

    % Initialize output variables
    t_iter_out = floor(t_iter/dt_out);
    phi_t = zeros(nx,ny,t_iter_out); 
    mass_t = zeros(t_iter_out,1);
    E_t = zeros(t_iter_out,1);
    
    phi_old = phi0; % Initialize chemical state
    r_old = r0_fun(phi0,hx,hy,C0); % Initialize sav state
    mass = ch_mass(phi_old,h2,nx,ny); % Initialize mass
    E = ch_discrete_energy(phi_old,h2,nx,ny,epsilon2); % Initialize energy
    
    % Store initial states into output variables
    [phi_t,mass_t,E_t] = store_data(phi_t,mass_t,E_t,phi_old,mass,E,1); 
    
%% Run SAV solver

    for i = 1:t_iter
        % Calculate current phi, r, mass and E
        [phi_new,r_new] = sav_solver(phi_old,r_old, ...
            hx,hy,k2,k4,dt,epsilon2, ...
            boundary,C0,Beta,gamma0);
        mass = ch_mass(phi_new,h2,nx,ny);
        E = ch_discrete_energy(phi_new,h2,nx,ny,epsilon2);

        % Store current states every dt_out time steps
        if mod(i,dt_out) == 0
            [phi_t,mass_t,E_t] = store_data(phi_t,mass_t,E_t,phi_new,mass,E,i/dt_out+1);
        end

        % Update iteration variables
        phi_old = phi_new;
        r_old = r_new;

        % Print percentage of completion
        if mod(i/t_iter*100,5) == 0
            fprintf('%3.0f percent complete\n',i/t_iter*100)
        end
    end

%% For post-processing

    if pathname == "cd"
        pathname = pwd;
    end
    % Print phi to file if specified
    if printphi
        Filename = strcat(pathname, 'phi.csv');
        % Path = strcat(pwd, '/', Filename);
        writematrix(phi_new, Filename, 'WriteMode', 'append'); 
        fprintf('Data appended to %s\n', Filename);
    end

    % Center mass and normalize energy to t == 0
    delta_mass_t = mass_t - mass_t(1);
    E_t = E_t/E_t(1);
    
    % Output t_out vector for post-processing
    t_out = (0:dt_out:t_iter)*dt;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions

    % Local function for calculating the discrete energy across the domain
    function E = ch_discrete_energy(phi,hxy,gridx,gridy,eps2)
        a = hxy*sum(sum(f(phi))); %Calculate chemical free energy
        sum_i = 0; % Initialize interfacial free energy in x
        for i = 1:gridx-1
            for j = 1:gridy
                sum_i = sum_i + (phi(i+1,j)-phi(i,j))^2;
            end
        end
        sum_j = 0; % Initialize interfacial free energy in y
        for i = 1:gridx
            for j = 1:gridy-1
                sum_j = sum_j + (phi(i,j+1)-phi(i,j))^2;
            end
        end
        E = a + 0.5*eps2*(sum_i+sum_j);
    end

    % Local function for calculating mass across the domain
    function mass = ch_mass(phi,h2,nx,ny)
        mass = sum(sum(phi))/(h2*nx*ny);
    end

    % Local function for data storage
    function [phi_t,mass_t,E_t] = store_data(phi_t,mass_t,E_t,phi,mass,E,step_i)
        phi_t(:,:,step_i) = phi;
        mass_t(step_i) = mass;
        E_t(step_i) = E;
    end


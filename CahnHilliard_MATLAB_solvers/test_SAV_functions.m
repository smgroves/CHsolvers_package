
% Make x matrix 2 x 2
x = [1.0, 2.0; 3.0, 4.0];
boundary = "neumann";
[nx,ny] = size(x); % Define number of grid points in x and y
xright = 1;
xleft = 0;
yright = 1;
yleft = 0;
Lx = xright-xleft; Ly = yright-yleft;

% Decide on the solver's mesh spacing for NEUMANN vs PERIODIC
%  - For Neumann: we will mirror the domain, so pass 2*hx and 2*hy into sav_solver.
%  - For Periodic: keep as-is.
if strcmpi(boundary,'neumann')
    Lx = 2*Lx;
    Ly = 2*Ly;
    nx = 2*nx;
    ny = 2*ny;
elseif strcmpi(boundary,'periodic')
    Lx = Lx;
    Ly = Ly;
    nx = nx;
    ny = ny;
end

hx = Lx/nx; hy = Ly/ny;
h2 = hx*hy; % Define mesh size

m = 8;
epsilon2 = NaN;
if isnan(epsilon2)
    epsilon2 = h2*m^2/(2*sqrt(2)*atanh(0.9))^2; % Define ϵ^2 if not prespecified
else
    m = sqrt((epsilon2*(2*sqrt(2)*atanh(0.9))^2)/h2); % Else overwrite m
    display(m);
end

k_x = 1i*[0:nx/2 -nx/2+1:-1]*(2*pi/Lx); k_y = 1i*[0:ny/2 -ny/2+1:-1]*(2*pi/Ly);
k_xx = k_x.^2; k_yy = k_y.^2;
[kxx,kyy] = meshgrid(k_xx,k_yy);
% k2 = kxx + kyy
% k4 = k2.^2



% r0_fun(x, hx, hy, 0)
% xext = ext(x)


% x_back = extback(xext)

% df(x,0)

% f(x)
% k2 .* fft2_filtered(xext)


% %
dt = 1e-5
gamma0 = 0;
Beta = 0;
C0 = 0;
r0 = r0_fun(ext(x),hx,hy,C0) % Initialize sav state

phi0_df   = df(ext(x),gamma0)        % df at phi0
Lap_dfphi0 = Lap_SAV(phi0_df, k2)                    % Lap of df(phi0)
phi_bar = A_inv_CN(ext(x) + dt/2 * Lap_dfphi0, dt, k2, k4, gamma0, epsilon2)

% Step 1
b = b_fun(phi_bar,hx,hy,C0,gamma0)
g = g_fun_CN(ext(x),r0,b,dt,hx,hy,epsilon2,gamma0,Beta,C0,k2)


function x_ext = ext(x)
    % Mirroring extension: takes x(nx, ny) -> x_ext(2*nx, 2*ny)
        [nx, ny] = size(x);
        x_ext = zeros(2*nx, 2*ny);
    
        % Original block
        x_ext(1:nx, 1:ny) = x;
    
        % Flip horizontally
        x_ext(1:nx, ny+1:2*ny) = x(:, end:-1:1);
    
        % Flip vertically
        x_ext(nx+1:2*nx, 1:ny) = x(end:-1:1, :);
    
        % Flip both
        x_ext(nx+1:2*nx, ny+1:2*ny) = x(end:-1:1, end:-1:1);
    end

% Local function for "flip" extension back
function x_back = extback(x_ext)
    % Shrinks from 2*nx x 2*ny back to nx x ny (upper-left block)
        [nx_ext, ny_ext] = size(x_ext);
        nx = nx_ext/2;
        ny = ny_ext/2;
        x_back = x_ext(1:nx, 1:ny);
end
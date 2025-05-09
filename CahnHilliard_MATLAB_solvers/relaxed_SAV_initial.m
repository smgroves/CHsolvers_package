% Define parameters
Lx = 1;          % Length in x-direction
Ly = 1;          % Length in y-direction
Nx = 128;        % Number of points in x-direction
Ny = 128;        % Number of points in y-direction
epsilon = 0.01;  % Small positive parameter (example value, can be adjusted)

% Create 1D vectors for x and y
x = linspace(0, 1, Nx);
y = linspace(0, 1, Ny);

% Create 2D grid
[X, Y] = meshgrid(x, y);

% Compute radial distance r
R = sqrt((X - Lx/2).^2 + (Y - Ly/2).^2);

% Compute angle theta
Theta = atan2(Y - 0.5*Ly, X - 0.5*Lx);

% Compute the function phi0
numerator = 1.5 + 1.2 * cos(6 * Theta) - 2 * pi * R;
phi0 = tanh(numerator / (sqrt(2) * epsilon));
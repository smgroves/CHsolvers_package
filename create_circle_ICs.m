% ----------------------------------------------------------
% create_circle_ICs.m  –  generate tanh-smoothed circles
% ----------------------------------------------------------
%  • Domain: [0,1] × [0,1]
%  • Radius : R = 0.1 (fraction of domain length)
%  • Interface thickness controlled by eps_fac
%  • Output : CSVs named  circle_R01_N{Nx}.csv  in ./IC
% ----------------------------------------------------------

clear; clc;

outDir  = 'IC';                          % save folder
if ~exist(outDir, 'dir'), mkdir(outDir); end

Nlist   = [64, 128, 256, 512];           % grid sizes
R       = 0.1;                           % physical radius
eps_fac = 0.01;                          % controls tanh width

for N = Nlist
    dx   = 1/N;                          % grid spacing
    eps  = eps_fac;                      % interface parameter
    [x,y] = meshgrid(dx/2:dx:1-dx/2);    % cell-centered grid
    r     = sqrt( (x-0.5).^2 + (y-0.5).^2 );

    % tanh profile: +1 inside, −1 outside
    phi = tanh( (R - r) / (sqrt(2)*eps) );

    % write to CSV
    filename = sprintf('circle_R01_N%d.csv', N);
    writematrix(phi, fullfile(outDir, filename));
    fprintf('Wrote %s\n', filename);
end

disp('All circle ICs generated.');
%RUN_FD_IMPLICIT_REVIEWER_TESTS Reviewer-facing implicit Euler FD controls.
%
% This script tests whether the finite-difference Cahn-Hilliard comparison
% remains stable when the legacy forward-Euler time step is replaced by a
% fully implicit backward-Euler step with the same spatial discretization.
%
% Outputs are written to:
%   output/implicit_fd_reviewer/

clear; clc; close all;

outdir = fullfile(pwd, 'output', 'implicit_fd_reviewer');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

cases = [
    struct('label', 'N64_full_gate', 'GridSize', 64,  'n_big', 200, 'do_small_consistency', true,  'seed', 1234)
    struct('label', 'N128_early',    'GridSize', 128, 'n_big', 20,  'do_small_consistency', false, 'seed', 1234)
];

summary = table();
results = struct();

for i = 1:numel(cases)
    fprintf('\n=== Running %s ===\n', cases(i).label);
    result = run_case(cases(i));
    results.(cases(i).label) = result; %#ok<SAGROW>
    if i == 1
        summary = result.summary;
    else
        summary = [summary; result.summary]; %#ok<AGROW>
    end
    plot_field_comparison(result, outdir);
    plot_diagnostics(result, outdir);
end

writetable(summary, fullfile(outdir, 'implicit_fd_reviewer_summary.csv'));
save(fullfile(outdir, 'implicit_fd_reviewer_results.mat'), 'results', 'summary');
disp(summary);

assert(all(summary.explicit_big_failed), ...
    'Large-step explicit FD did not fail or become physically unusable in one or more test cases.');
assert(all(summary.implicit_big_newton_all_ok), ...
    'Implicit Euler Newton iteration failed in one or more test cases.');
assert(max(summary.implicit_big_newton_res_max) < 1e-8, ...
    'Implicit Euler maximum Newton residual is larger than expected.');

assert(summary.small_dt_explicit_vs_implicit_RMSE(1) < 5e-3, ...
    '64^2 small-dt implicit/explicit RMSE is unexpectedly large.');
assert(summary.implicit_big_mass_error(1) < 1e-10, ...
    '64^2 implicit Euler did not conserve mass to tolerance.');
assert(summary.implicit_big_energy_ratio(1) <= 1.0 + 1e-10, ...
    '64^2 implicit Euler energy did not decrease.');
assert(summary.implicit_big_vs_explicit_reference_RMSE(1) < 5e-3, ...
    '64^2 implicit large-dt solution differs too much from the explicit reference.');
assert(summary.implicit_big_vs_explicit_reference_corr(1) > 0.999, ...
    '64^2 implicit large-dt solution correlation with reference is too low.');

assert(summary.implicit_big_mass_error(2) < 1e-10, ...
    '128^2 implicit Euler did not conserve mass to tolerance.');
assert(summary.implicit_big_energy_ratio(2) <= 1.0 + 1e-10, ...
    '128^2 implicit Euler energy did not decrease.');
assert(summary.implicit_big_vs_explicit_reference_RMSE(2) < 5e-3, ...
    '128^2 implicit early-time solution differs too much from the explicit reference.');
assert(summary.implicit_big_vs_explicit_reference_corr(2) > 0.999, ...
    '128^2 implicit early-time solution correlation with reference is too low.');

fprintf('\nReviewer implicit-FD tests passed. Results written to %s\n', outdir);

function result = run_case(cfg)
boundary = 'neumann';
m = 8;
gamma = (m / (2 * sqrt(2) * atanh(0.9)))^2;
D = cfg.GridSize^2;
dt_big = 5.5e-6 * (128 / cfg.GridSize)^2;
dt_small = dt_big / 100;
n_big = cfg.n_big;
n_ref = n_big * 100;

phi0 = ch_make_spinodal_ic(cfg.GridSize, 0.50, ...
    'seed', cfg.seed, 'boundary', boundary, 'n_relax', 4, 'smooth_theta', 0.1);

small_rmse = NaN;
small_mass_err = NaN;
u_exp_small = NaN(size(phi0));
u_imp_small = NaN(size(phi0));

if cfg.do_small_consistency
    n_consistency = 100;
    [phi_exp_small, st_exp_small] = ch_fd_solve(phi0, ...
        'method', 'explicit', 'dt', dt_small, 't_iter', n_consistency, ...
        'D', D, 'gamma', gamma, 'boundary', boundary, 'save_every', n_consistency);
    [phi_imp_small, st_imp_small] = ch_fd_solve(phi0, ...
        'method', 'implicit_euler', 'dt', dt_small, 't_iter', n_consistency, ...
        'D', D, 'gamma', gamma, 'boundary', boundary, 'save_every', n_consistency, ...
        'newton_tol', 1e-10, 'newton_max_iter', 20, 'linear_solver', 'backslash');
    u_exp_small = phi_exp_small(:, :, end);
    u_imp_small = phi_imp_small(:, :, end);
    small_rmse = rmse(u_exp_small, u_imp_small);
    small_mass_err = max(abs(st_imp_small.mass - st_imp_small.mass(1)));
else
    st_exp_small = [];
    st_imp_small = [];
end

[phi_exp_big, st_exp_big] = ch_fd_solve(phi0, ...
    'method', 'explicit', 'dt', dt_big, 't_iter', n_big, ...
    'D', D, 'gamma', gamma, 'boundary', boundary, 'save_every', 1);
explicit_big_last = phi_exp_big(:, :, end);
explicit_big_max_abs = max(abs(explicit_big_last(:)));

[phi_imp_big, st_imp_big] = ch_fd_solve(phi0, ...
    'method', 'implicit_euler', 'dt', dt_big, 't_iter', n_big, ...
    'D', D, 'gamma', gamma, 'boundary', boundary, 'save_every', n_big, ...
    'newton_tol', 1e-9, 'newton_max_iter', 20, 'linear_solver', 'backslash');
implicit_big_final = phi_imp_big(:, :, end);
implicit_mass_err = max(abs(st_imp_big.mass - st_imp_big.mass(1)));
implicit_energy_ratio = st_imp_big.energy(end) / st_imp_big.energy(1);
implicit_newton_mean = mean(st_imp_big.newton_iter, 'omitnan');
implicit_newton_max = max(st_imp_big.newton_iter);
implicit_newton_all_ok = all(st_imp_big.newton_ok);
implicit_newton_res_max = max(st_imp_big.newton_res);
explicit_big_failed = ~isnan(st_exp_big.blowup_step) || explicit_big_max_abs > 10;

[phi_exp_ref, st_exp_ref] = ch_fd_solve(phi0, ...
    'method', 'explicit', 'dt', dt_small, 't_iter', n_ref, ...
    'D', D, 'gamma', gamma, 'boundary', boundary, 'save_every', n_ref);
reference_final = phi_exp_ref(:, :, end);
big_vs_ref_rmse = rmse(implicit_big_final, reference_final);
cc = corrcoef(implicit_big_final(:), reference_final(:));
big_vs_ref_corr = cc(1, 2);

result.label = cfg.label;
result.GridSize = cfg.GridSize;
result.dt_big = dt_big;
result.dt_small = dt_small;
result.n_big = n_big;
result.n_ref = n_ref;
result.phi0 = phi0;
result.explicit_small_final = u_exp_small;
result.implicit_small_final = u_imp_small;
result.explicit_big_last = explicit_big_last;
result.implicit_big_final = implicit_big_final;
result.reference_final = reference_final;
result.st_exp_small = st_exp_small;
result.st_imp_small = st_imp_small;
result.st_exp_big = st_exp_big;
result.st_imp_big = st_imp_big;
result.st_exp_ref = st_exp_ref;

result.summary = table( ...
    string(cfg.label), cfg.GridSize, dt_big, dt_small, n_big, n_ref, ...
    small_rmse, small_mass_err, st_exp_big.blowup_step, explicit_big_max_abs, explicit_big_failed, ...
    implicit_mass_err, implicit_energy_ratio, implicit_newton_mean, implicit_newton_max, ...
    implicit_newton_all_ok, implicit_newton_res_max, ...
    big_vs_ref_rmse, big_vs_ref_corr, st_imp_big.elapsed, st_exp_ref.elapsed, ...
    'VariableNames', {'case_label', 'GridSize', 'dt_big', 'dt_small', 'n_big', 'n_ref', ...
    'small_dt_explicit_vs_implicit_RMSE', 'small_dt_implicit_mass_error', ...
    'explicit_big_blowup_step_NaN_if_not_stopped', 'explicit_big_max_abs_phi', ...
    'explicit_big_failed', ...
    'implicit_big_mass_error', 'implicit_big_energy_ratio', ...
    'implicit_big_newton_iter_mean', 'implicit_big_newton_iter_max', ...
    'implicit_big_newton_all_ok', 'implicit_big_newton_res_max', ...
    'implicit_big_vs_explicit_reference_RMSE', 'implicit_big_vs_explicit_reference_corr', ...
    'implicit_big_elapsed_seconds', 'explicit_reference_elapsed_seconds'});
end

function plot_field_comparison(result, outdir)
fig = figure('visible', 'off', 'Position', [100 100 1500 560]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_field_tile(result.phi0, 'Initial condition', [-1 1]);
plot_field_tile(result.reference_final, 'Explicit reference', [-1 1]);
plot_field_tile(result.implicit_big_final, 'Implicit Euler FD', [-1 1]);

comparison_rmse = rmse(result.implicit_big_final, result.reference_final);
cc = corrcoef(result.implicit_big_final(:), result.reference_final(:));
comparison_corr = cc(1, 2);

colormap(fig, redblue(256));
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = '\phi';
cb.FontSize = 12;
title(tl, sprintf('%s final-field agreement: RMSE = %.3g, corr = %.6f', ...
    result.label, comparison_rmse, comparison_corr), ...
    'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold');
exportgraphics(fig, fullfile(outdir, sprintf('%s_fields.png', result.label)), 'Resolution', 200);
close(fig);
end

function plot_diagnostics(result, outdir)
fig = figure('visible', 'off', 'Position', [100 100 1000 620]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(result.st_imp_big.t, result.st_imp_big.energy_norm, 'LineWidth', 1.7); hold on;
plot(result.st_exp_ref.t, result.st_exp_ref.energy_norm, '--', 'LineWidth', 1.4);
xlabel('time');
ylabel('E/E_0');
ylim([0, 1.05]);
legend({'implicit large dt', 'explicit tiny-dt reference'}, 'Location', 'best');
title(sprintf('%s energy; explicit large dt fails at step %.0f', ...
    result.label, result.st_exp_big.blowup_step), 'Interpreter', 'none');
grid on;

nexttile;
plot(result.st_imp_big.t, result.st_imp_big.mass - result.st_imp_big.mass(1), 'LineWidth', 1.7); hold on;
plot(result.st_exp_ref.t, result.st_exp_ref.mass - result.st_exp_ref.mass(1), '--', 'LineWidth', 1.4);
xlabel('time');
ylabel('mass error');
legend({'implicit large dt', 'explicit tiny-dt reference'}, 'Location', 'best');
title('mass conservation');
grid on;

exportgraphics(fig, fullfile(outdir, sprintf('%s_diagnostics.png', result.label)), 'Resolution', 200);
close(fig);
end

function plot_field_tile(phi, ttl, lims)
nexttile;
imagesc(phi);
axis image off;
title(ttl, 'Interpreter', 'none', 'FontSize', 15, 'FontWeight', 'bold');
clim(lims);
end

function value = rmse(a, b)
d = a(:) - b(:);
value = sqrt(mean(d.^2));
end

function c = redblue(m)
if nargin < 1
    m = size(get(gcf, 'colormap'), 1);
end

if mod(m, 2) == 0
    m1 = m * 0.5;
    r = (0:m1-1)' / max(m1-1, 1);
    g = r;
    r = [r; ones(m1, 1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    m1 = floor(m * 0.5);
    r = (0:m1-1)' / max(m1, 1);
    g = r;
    r = [r; ones(m1+1, 1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b];
end

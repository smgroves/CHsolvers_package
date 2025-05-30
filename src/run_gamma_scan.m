function out = run_gamma_scan(phi0, dtList, gammaList, ...
                              Tfinal, boundary, eta)
% Sweep over stabilisation parameter γ₀ (eta fixed).
% Returns struct array: out(k).gamma0 , out(k).table

nG  = numel(gammaList);
out = struct([]);

for k = 1:nG
    g0 = gammaList(k);
    fprintf('--- γ₀ = %.3g ---\n', g0);

    param.gamma0   = g0;
    param.eta      = eta;
    param.boundary = boundary;

    out(k).gamma0 = g0;
    out(k).table  = benchmark_SAV(phi0, dtList, Tfinal, param);
end
end
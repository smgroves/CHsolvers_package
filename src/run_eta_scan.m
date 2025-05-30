function out = run_eta_scan(phi0, dtList, etaList, ...
                            Tfinal, boundary, gamma0)
% Sweep over relaxation parameter η.
% Uses benchmark_SAV(phi0, dtList, Tfinal, param)
% Returns struct array: out(k).eta   , out(k).table

nEta = numel(etaList);
out  = struct([]);

for k = 1:nEta
    eta = etaList(k);
    fprintf('--- η = %.3g ---\n', eta);

    param.gamma0   = gamma0;
    param.eta      = eta;
    param.boundary = boundary;

    out(k).eta   = eta;
    out(k).table = benchmark_SAV(phi0, dtList, Tfinal, param);
end
end
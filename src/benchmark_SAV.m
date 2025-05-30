function T = benchmark_SAV(phi0,dtList,Tfinal,param)
% Run CahnHilliard_SAV at each Δt until t = Tfinal.
% Reference solution = finest Δt.
% Returns a table with columns: dt, L2

nDt      = numel(dtList);
phiRef   = [];
dtRec    = zeros(nDt,1);
l2Rec    = zeros(nDt,1);

% --- reference run (smallest Δt) ---
dt_ref  = dtList(end);
iterRef = ceil(Tfinal/dt_ref);
dt_out  = max(1,round(iterRef/200));          % ~200 outputs

[~,phi_t,~,~,~] = CahnHilliard_SAV( ...
    phi0,'t_iter',iterRef,'dt',dt_ref,'dt_out',dt_out, ...
    'boundary',param.boundary,'gamma0',param.gamma0, ...
    'eta',param.eta,'m',8,'xi_flag',1,'C0',0);
phiRef = phi_t(:,:,end);

% --- coarser Δt loop ---
for j = 1:nDt
    dt      = dtList(j);
    nSteps  = ceil(Tfinal/dt);
    dt_out  = max(1,round(nSteps/200));

    fprintf('  Δt = %.1e --> %d steps ... ',dt,nSteps);

    [~,phi_t,~,~,~] = CahnHilliard_SAV( ...
        phi0,'t_iter',nSteps,'dt',dt,'dt_out',dt_out, ...
        'boundary',param.boundary,'gamma0',param.gamma0, ...
        'eta',param.eta,'m',8,'xi_flag',1,'C0',0);

    l2Rec(j) = compute_L2(phi_t(:,:,end),phiRef);
    dtRec(j) = dt;

    fprintf('L2 = %.3e\n',l2Rec(j));
end

T = table(dtRec,l2Rec,'VariableNames',{'dt','L2'});
end
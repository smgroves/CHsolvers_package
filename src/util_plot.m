function util_plot(resultStruct,modeFlag)
% Two–panel figure: (a)  L2 vs Δt
figure('Name',['Convergence-' modeFlag],'Position',[120 120 650 500]);

subplot(1,1,1); hold on; grid on;
xlabel('\Delta t'); ylabel('L_2 error');
set(gca,'XScale','log','YScale','log');

leg = {};
for k = 1:numel(resultStruct)
    T = resultStruct(k).table;
    plot(T.dt,T.L2,'-o','LineWidth',1.4);
    if isfield(resultStruct(k),'eta')
        leg{end+1} = sprintf('\\eta = %.3g',resultStruct(k).eta);
    else
        leg{end+1} = sprintf('\\gamma_0 = %.3g',resultStruct(k).gamma0);
    end
end
legend(leg,'Location','SouthWest');
title(['SAV convergence – ' modeFlag]);
end
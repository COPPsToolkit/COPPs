cd('Models\SW07')

%% Generate SW07 IRFs

% Standard Phillips curve
dynare Smets_Wouters_2007_GB09 -DIRFs noclearall
load Smets_Wouters_2007_GB09_results
oo_SWorg = oo_;
save oo_SWorg_file oo_SWorg

% Flat Phillips curve
dynare Smets_Wouters_2007_GB09 -DIRFs -DFlatPC
load Smets_Wouters_2007_GB09_results
load oo_SWorg_file

%% Read FRB/US IRFs

[XlsValues XlsText] = xlsread('FRBUS_IRFs_2018','FFR_VAR');

% save VAR-based expectations
Frbus.VAR.previous = XlsValues(:,[1 3 7]);
Frbus.VAR.new = XlsValues(:,[2 4 8]);

[XlsValues XlsText] = xlsread('FRBUS_IRFs_2018','FFR_MCE');

% save Model Consistent Expectations (MCE)
Frbus.MCE.previous = XlsValues(:,[1 3 7]);
Frbus.MCE.new = XlsValues(:,[2 4 8]);

%%

Adjustment = 0.97 / max(oo_.irfs.obs_r_ann_em);

figure
subplot(2,2,2)
hold all; grid on
plot(Frbus.VAR.new(1:40,1), 'g-','LineWidth',2)
plot(Frbus.MCE.new(1:40,1), 'g-.','LineWidth',2)
plot(Frbus.VAR.previous(1:40,1), 'k-','LineWidth',2)
plot(Frbus.MCE.previous(1:40,1), 'k-.','LineWidth',2)
plot(oo_SWorg.irfs.ygap_em * Adjustment, 'r--','LineWidth',2)
plot(oo_.irfs.ygap_em * Adjustment, 'b--','LineWidth',2)
title('Output gap')
ylim([-.6 0.2])
subplot(2,2,1)
hold all; grid on
plot(Frbus.VAR.new(1:40,2), 'g-','LineWidth',2)
plot(Frbus.MCE.new(1:40,2), 'g-.','LineWidth',2)
plot(Frbus.VAR.previous(1:40,2), 'k-','LineWidth',2)
plot(Frbus.MCE.previous(1:40,2), 'k-.','LineWidth',2)
plot(oo_SWorg.irfs.pinf_4q_em * Adjustment, 'r--','LineWidth',2)
plot(oo_.irfs.pinf_4q_em * Adjustment, 'b--','LineWidth',2)
title('Inflation')
ylim([-0.25 0.05])
subplot(2,2,3)
hold all; grid on
plot(Frbus.VAR.new(1:40,3), 'g-','LineWidth',2)
plot(Frbus.MCE.new(1:40,3), 'g-.','LineWidth',2)
plot(Frbus.VAR.previous(1:40,3), 'k-','LineWidth',2)
plot(Frbus.MCE.previous(1:40,3), 'k-.','LineWidth',2)
plot(oo_SWorg.irfs.obs_r_ann_em * Adjustment, 'r--','LineWidth',2)
plot(oo_.irfs.obs_r_ann_em * Adjustment, 'b--','LineWidth',2)
title('Interest rate')
ylim([-0.5 1.0])
legend('New FRBUS - VAR expectations','New FRBUS - MCE','Old FRBUS - VAR expectations','Old FRBUS - MCE','SW07: Steep PC','SW07: Flat PC')

%%
cd ..
cd ..
cd ..
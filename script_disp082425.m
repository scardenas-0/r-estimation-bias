path(path,'..')

%%
% Figure for showing what combinations of parameters would be compatible
% with observations of measles

d1 = load('Data/082225_msls_r.mat');
% d2 = load('Data/061925_msls_02_perf');
R_inf = squeeze(d1.rs_inf_res)';

figure(1),clf

%Third colm of d2.res_arr: [C_pp C_ps C_sp, C_ss perf_p perf_s accur r_inf r_bias]

% From McQuillan et al., The Journal of Infectious Diseases, Volume 196, Issue 10, 15 November 2007, Pages 1459?1464
% During 1999-2004, the rate of measles seropositivity in the population overall was 95.9% (95% confidence interval [CI], 95.1%?96.5%)
%
% Based on this, we'll assume the baseline R is based on a uniform
% seropositivity of 95.9
% **updated to be 95%
%
% R0 = R_inf/.041;
% Assumign 100% of vaccination leads to seroprevalene, critical decrease in vaccination satisfies
% 1 = (.05+x)*R0 = (1 + x/.041)*R_inf
% Thus x = .05(1/R_inf -1)*100

% contour(d1.rp_arr,d1.p_obs_arr,4.1*(1./R_inf -1),[1.1, 1.2, 1.4, 1.6, 1.8, 2],'ShowText','On','LineWidth',2)
contour(d1.rp_arr,d1.p_obs_arr,5.562864779*(1./R_inf -1),'ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Absolute decrease in vaccination percentage','leading to critical transmission threshold'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'Primary cases per cluster (R_p)'})
%tt = text(2.01, .9, {'Heterogenous','R_s = 0.8'})
%set(tt, 'FontSize', 10, 'FontWeight','bold')
%ylim([.5 1])

% savepdf('Figures/040519_msls_app')
exportgraphics(gcf, 'Figures/082225_msls_vax_app.pdf')


%%
% Figure look at R and OR for measles in 2020s
figure(1),clf

clear
d2 = load('Data/082225_msls_OR.mat');
d1 = d2.d_msls;
R_inf = squeeze(d1.rs_inf_res)';

% Prob observed primary is a primary
subplot (3,2,1)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,5)',[0.7, 0.8, 0.9, 1],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a case classified_ ','as primary is primary (P_{p\rightarrow p})'})
%title({'P_{p\rightarrow p}'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})

% Prob observed secondary is a secondary
subplot (3,2,3)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,6)',[0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a case classified_ ','as secondary is secondary (P_{s\rightarrow s})'})
%title({'P_{s\rightarrow s}'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})
ylabel({'Observation probability (P_{obs})'},'FontSize',11)

% Inferred proportion positive among primary
% res_arr(:,:,10:13) = [Q_pn Q_pp Q_sn Q_sp];
subplot (3,2,2)
temp = d2.res_arr(:,:,11)./(d2.res_arr(:,:,10) + d2.res_arr(:,:,11));
contour(d1.rp_arr,d1.p_obs_arr,temp',[0.22, 0.24, 0.26, 0.28, 0.3],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Inferred proportion of','primary cases that are >19 years old'})
%xlabel({'              Primary cases per cluster (R_p)'})

% Inferred proportion positive among secondary
% res_arr(:,:,10:13) = [Q_pn Q_pp Q_sn Q_sp];
subplot (3,2,4)
temp = d2.res_arr(:,:,13)./(d2.res_arr(:,:,12) + d2.res_arr(:,:,13));
contour(d1.rp_arr,d1.p_obs_arr,temp','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Inferred proportion of','secondary cases that are >19 years old'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})

% Inferred r_effective
subplot (3,2,5)
contour(d1.rp_arr,d1.p_obs_arr,R_inf,[0.55, 0.6, 0.65, 0.7, 0.75, 0.8],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Effective reproduction number'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel('Primaries per cluster (R_p)')
text(1.0,-0.08,0,'Primary cases per cluster (R_p)','FontSize', 11)

% Inferred odds ratio
% res_arr(:,:,10:13) = [Q_pn Q_pp Q_sn Q_sp];
subplot (3,2,6)
temp = d2.res_arr(:,:,11).*d2.res_arr(:,:,12)./(d2.res_arr(:,:,10) .* d2.res_arr(:,:,13));
% bounds(temp)
% abs(temp) used because of some convergence difficulties yielding large negative
% values of OR for high p_obs and r_p
% contour(d1.rp_arr,d1.p_obs_arr,abs(temp)',[0.01 0.02 0.03 0.04 0.05 0.06],'ShowText','On','LineWidth',2)
contour(d1.rp_arr,d1.p_obs_arr,temp',[4, 4.5, 5, 5.5, 6],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Inferred odds ratio for','>19 year old cases being primary'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel('Primaries per cluster (R_p)')
text(1.0,-0.08,0,'Primary cases per cluster (R_p)','FontSize', 11)

% savepdf('Figures/040519_mpx_app')
set(gcf, 'position', [10, 10, 550, 600])
exportgraphics(gcf, 'Figures/082425_msls_app.pdf')


%% Figure look at R and OR for mpox in 2020s
figure(1),clf

clear
d2 = load('Data/082225_mpox_OR.mat');
d1 = d2.d_mpox;
R_inf = squeeze(d1.rs_inf_res)';

% Prob observed primary is a primary
subplot (3,2,1)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,5)',[0.7, 0.8, 0.9, 1],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a case classified_ ','as primary is primary (P_{p\rightarrow p})'})
%title({'P_{p\rightarrow p}'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})

% Prob observed secondary is a secondary
subplot (3,2,3)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,6)','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a case classified_ ','as secondary is secondary (P_{s\rightarrow s})'})
%title({'P_{s\rightarrow s}'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})
ylabel({'Observation probability (P_{obs})'},'FontSize',11)

% Inferred proportion positive among primary
% res_arr(:,:,10:13) = [Q_pn Q_pp Q_sn Q_sp];
subplot (3,2,2)
temp = d2.res_arr(:,:,11)./(d2.res_arr(:,:,10) + d2.res_arr(:,:,11));
contour(d1.rp_arr,d1.p_obs_arr,temp',[0.55, 0.56, 0.57, 0.58],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Inferred proportion of primary','cases that report animal contact'})

%xlabel({'              Primary cases per cluster (R_p)'})

% Inferred proportion positive among secondary
% res_arr(:,:,10:13) = [Q_pn Q_pp Q_sn Q_sp];
subplot (3,2,4)
temp = d2.res_arr(:,:,13)./(d2.res_arr(:,:,12) + d2.res_arr(:,:,13));
contour(d1.rp_arr,d1.p_obs_arr,temp',[0.42, 0.43, 0.44, 0.45],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Inferred proportion of secondary','cases that report animal contact'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})

% Inferred r_effective
subplot (3,2,5)
contour(d1.rp_arr,d1.p_obs_arr,R_inf,'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Effective reproduction number'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel('Primaries per cluster (R_p)')
text(1.0,-0.08,0,'Primary cases per cluster (R_p)','FontSize', 11)

% Inferred odds ratio
% res_arr(:,:,10:13) = [Q_pn Q_pp Q_sn Q_sp];
subplot (3,2,6)
temp = d2.res_arr(:,:,11).*d2.res_arr(:,:,12)./(d2.res_arr(:,:,10) .* d2.res_arr(:,:,13));
% bounds(temp)
% abs(temp) used because of some convergence difficulties yielding large negative
% values of OR for high p_obs and r_p
% contour(d1.rp_arr,d1.p_obs_arr,abs(temp)',[0.01 0.02 0.03 0.04 0.05 0.06],'ShowText','On','LineWidth',2)
contour(d1.rp_arr,d1.p_obs_arr,temp',[1.5, 1.55, 1.6, 1.65, 1.7],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Inferred odds ratio for primary','cases reporting animal contact'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel('Primaries per cluster (R_p)')
text(1.0,-0.08,0,'Primary cases per cluster (R_p)','FontSize', 11)

% savepdf('Figures/040519_mpx_app')
set(gcf, 'position', [10, 10, 550, 600])
exportgraphics(gcf, 'Figures/111925_mpox_app.pdf')

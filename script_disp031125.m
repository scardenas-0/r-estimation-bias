path(path,'..')
%%
% Figure for describing how often primary cases are classified as primary
% vs secondary

d1 = load('Data/011119_rp_pobs');

figure(1),clf
set(gcf, 'Position', [10, 10, 1200, 900])
subplot(2,2,1)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,1))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that a primary case is_ ', 'classified as a primary case (C_{p\rightarrow p})'})
%xlabel({'              Primary cases per cluster (R_p)'})
ylabel({'Observation_ ','probability (P_{obs})'},'FontSize',16)
subplot(2,2,2)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,2))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
%xlabel({'              Primary cases per cluster (R_p)'})
title({'Probability that a primary case is_ ', 'classified as a secondary case (C_{p\rightarrow s})'})
set(gca,'YTickLabel','')
text(.5,-.05,0,'Primary cases per cluster (R_p)','FontSize', 16)

% savepdf('Figures/class_prim040519')
% saveas(gcf,'Figures/class_prim012518','pdf')
exportgraphics(gcf, 'Figures/class_prim040519.pdf')

% rp,rs,k,p,[C_pp C_ps C_sp, C_ss perf_p perf_s accur r_inf r_bias]
%%
% Figure for describing how often secondary cases are classified as primary
% vs secondary

d1 = load('Data/011119_rp_pobs');
figure(1),clf
subplot(3,2,1)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,3))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that a secondary case is_ ', 'classified as a primary case (C_{s\rightarrow p})'})
%ylabel({'Observation_ ','probability (P_{obs})'})
set(gca,'XTickLabel','')
tt = text(2.01, .9, {'Heterogeneous','R_s = 0.2'})
set(tt, 'FontSize', 10, 'FontWeight','bold')
subplot(3,2,2)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,4))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that a secondary case is_ ', 'classified as a secondary case (C_{s\rightarrow s})'})
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')

subplot(3,2,3)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,3,1,:,3))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
ylabel({'Observation probability (P_{obs})'},'FontSize',16)
set(gca,'XTickLabel','')
tt = text(2.01, .9, {'Heterogeneous','R_s = 0.8'})
set(tt, 'FontSize', 10, 'FontWeight','bold')
subplot(3,2,4)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,3,1,:,4))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')

subplot(3,2,5)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,3,:,3))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})
tt = text(2.01, .9, {'Homogeneous','R_s = 0.2'})
set(tt, 'FontSize', 10, 'FontWeight','bold')
subplot(3,2,6)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,3,:,4))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
%xlabel({'              Primary cases per cluster (R_p)'})
set(gca,'YTickLabel','')
text(.5,-.1,0,'Primary cases per cluster (R_p)','FontSize', 16)

% savepdf('Figures/class_sec040519')
exportgraphics(gcf, 'Figures/class_sec040519.pdf')
% rp,rs,k,p,[C_pp C_ps C_sp, C_ss perf_p perf_s accur r_inf r_bias]

%%
% Figure for describing how trustworthy the results of null clasification
% are
d1 = load('Data/011119_rp_pobs');
figure(1),clf
subplot(3,3,1)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,5))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that a case_ ','classified as a primary case_ ', 'is a primary case (P_{p\rightarrow p})'})
%ylabel({'Observation_ ','probability (P_{obs})'})
set(gca,'XTickLabel','')
subplot(3,3,2)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,6))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that a case_ ','classified as a secondary case_ ', 'is a secondary case (P_{s\rightarrow s})'})
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
subplot(3,3,3)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,7))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that an observed_ ','case is classified correctly (\Theta)'})
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
tt = text(2.01, .9, {'Heterogeneous','R_s = 0.2'})
set(tt, 'FontSize', 10, 'FontWeight','bold')

subplot(3,3,4)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,3,1,:,5))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
ylabel({'Observation probability (P_{obs})'},'FontSize', 16)
set(gca,'XTickLabel','')
subplot(3,3,5)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,3,1,:,6))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
subplot(3,3,6)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,3,1,:,7))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
tt = text(2.01, .9, {'Heterogeneous','R_s = 0.8'})
set(tt, 'FontSize', 10, 'FontWeight','bold')

subplot(3,3,7)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,3,:,5))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
%xlabel({'Primary cases per cluster (R_p)'})
%ylabel({'Observation_ ','probability (P_{obs})'})
subplot(3,3,8)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,3,:,6))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
xlabel({'Primary cases per cluster (R_p)'},'FontSize', 16)
set(gca,'YTickLabel','')
subplot(3,3,9)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,3,:,7))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
%xlabel({'Primary cases per cluster (R_p)'})
set(gca,'YTickLabel','')
tt = text(2.01, .9, {'Homogeneous','R_s = 0.2'})
set(tt, 'FontSize', 10, 'FontWeight','bold')

% savepdf('Figures/class_accuracy040519')
exportgraphics(gcf, 'Figures/class_accuracy040519.pdf')
% rp,rs,k,p,[C_pp C_ps C_sp, C_ss perf_p perf_s accur r_inf r_bias]
%%
% Figure for describing how well classifier is able to estimate R

d1 = load('Data/011119_rp_pobs');
figure(1),clf
subplot(3,2,1)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,9))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Bias of R_s estimate (\delta_R)'})
%ylabel({'Observation_ ','probability (P_{obs})'})
set(gca,'XTickLabel','')
tt = text(2.01, .9, {'R_s = 0.2'})
set(tt, 'FontSize', 10, 'FontWeight','bold')
subplot(3,2,3)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,3,1,:,9))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
%ylabel({'Observation_ ','probability (P_{obs})'})
xlabel({'Primary cases per cluster (R_p)'},'FontSize', 16)
tt = text(2.01, .9, {'R_s = 0.8'})
set(tt, 'FontSize', 10, 'FontWeight','bold')
text(0.86,.5,0,'Primary cases per cluster (R_p)','FontSize', 16,'Rotation', 90)

% savepdf('Figures/R_est040519')
exportgraphics(gcf, 'Figures/R_est040519.pdf')
% rp,rs,k,p,[C_pp C_ps C_sp, C_ss perf_p perf_s accur r_inf r_bias]
%%
% Figure Bias of observed odds ratio
figure(1),clf
d = load('Data/040419_ORx_example')

nn = 0;
for rr = 1:length(d.rs_arr)
    for ii = 1:length(d.inc_arr)
        nn = nn+1;
        subplot(length(d.rs_arr),length(d.inc_arr),nn)
        contour(d.rp_arr,d.p_obs_arr,squeeze(d.res_arr(:,:,rr,ii,2)'), [1.5 2 2.5 3 3.5],'ShowText','On','LineWidth',2)
        title(strcat(['R_s =',' ',num2str(d.rs_arr(rr)),'   I_p =',' ',num2str(d.inc_arr(ii,1))]));
        set(gca,'FontSize',12)
        switch nn
            case 1
                set(gca,'XTickLabel','')
            case 2
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
            case 3
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
            case 4
                set(gca,'XTickLabel','')
                ylabel('Observaton probability (P_{obs})','FontSize',16)
            case 5
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
            case 6
                set(gca,'XTickLabel','')
                set(gca,'YTickLabel','')
            case 7
            case 8
                set(gca,'YTickLabel','')
                xlabel('Primaries per cluster (R_p)','FontSize',16)
            case 9
                set(gca,'YTickLabel','')
        end
    end
end

% savepdf('Figures/040519_ORx')
exportgraphics(gcf, 'Figures/040519_ORx.pdf')
%%
% Figure for showing what combinations of parameters would be compatible
% with observations of measles

%d1 = load('Data/010919_msls_model');
d1 = load('Data/043025_msls_model');
d2 = load('Data/043025_msls_perf');
R_inf = squeeze(d1.rs_inf_res)';

figure(1),clf

%Third colm of d2.res_arr: [C_pp C_ps C_sp, C_ss perf_p perf_s accur r_inf r_bias]


% Prob primary cases are labeled as primary
subplot (3,2,1)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,1)','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'C_{p\rightarrow p}'})
ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})


% Prob secondary cases are labeled as secondary
subplot (3,2,2)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,4)','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'C_{s\rightarrow s}'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})

clf
% Prob observed primary is a primary
%subplot (3,2,3)
subplot (2,2,1)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,5)','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'P_{p\rightarrow p}'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})

% Prob observed secondary is a secondary
%subplot (3,2,4)
subplot (2,2,2)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,6)','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'P_{s\rightarrow s}'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})


% From McQuillan et al., The Journal of Infectious Diseases, Volume 196, Issue 10, 15 November 2007, Pages 1459?1464
% During 1999-2004, the rate of measles seropositivity in the population overall was 95.9% (95% confidence interval [CI], 95.1%?96.5%)
%
% Based on this, we'll assume the baseline R is based on a uniform
% seropositivity of 95.9
% **updated to be 91.6%
%
% R0 = R_inf/.041;
% Assumign 100% of vaccination leads to seroprevalene, critical decrease in vaccination satisfies
% 1 = (.041+x)*R0 = (1 + x/.041)*R_inf
% Thus x = .041(1/R_inf -1)*100

% R effective
%subplot (3,2,5)
subplot (2,2,3)
contour(d1.rp_arr,d1.p_obs_arr,R_inf,'ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Effective reproduction number'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'Primary cases per cluster (R_p)'})
text(.96,.90,0,'Observation probability (P_{obs})','FontSize', 16,'Rotation', 90)
text(1.25,0.15,0,'Primary cases per cluster (R_p)','FontSize', 16)

%subplot (3,2,6)
subplot (2,2,4)
% contour(d1.rp_arr,d1.p_obs_arr,4.1*(1./R_inf -1),[1.1, 1.2, 1.4, 1.6, 1.8, 2],'ShowText','On','LineWidth',2)
contour(d1.rp_arr,d1.p_obs_arr,8.4*(1./R_inf -1),[3, 4, 5, 6, 7, 8, 9, 10, 11],'ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Absolute decrease in vaccination percentage','leading to critical transmission threshold'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'Primary cases per cluster (R_p)'})
%tt = text(2.01, .9, {'Heterogenous','R_s = 0.8'})
%set(tt, 'FontSize', 10, 'FontWeight','bold')
%ylim([.5 1])

% subplot (2,3,4)
% contour(d1.rp_arr,d1.p_obs_arr,squeeze(compat_arr)','ShowText','On','LineWidth',2)
% set(gca,'FontSize',12)
% title({'Bias of R_s estimate for measles'})
% ylabel({'Observation_ ','probability (P_{obs})'})
% xlabel({'              Primary cases per cluster (R_p)'})
% %tt = text(2.01, .9, {'Heterogenous','R_s = 0.8'})
% %set(tt, 'FontSize', 10, 'FontWeight','bold')
% ylim([.7 1])

% savepdf('Figures/040519_msls_app')
exportgraphics(gcf, 'Figures/052125_msls_34_app.pdf')

%%
% Figure look at R and OR for MPX in 80s
figure(1),clf

clear
d2 = load('Data/050125_nipah_OR');
% d1 = d2.d_mpx
d1 = d2.d_nipah
R_inf = squeeze(d1.rs_inf_res)';

% Prob observed primary is a primary
subplot (3,2,1)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,5)','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that a case_ ','classified as a primary case_ ', 'is a primary case (P_{p\rightarrow p})'})
%title({'P_{p\rightarrow p}'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})

% Prob observed secondary is a secondary
subplot (3,2,2)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,6)','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that a case_ ','classified as a secondary case_ ', 'is a secondary case (P_{s\rightarrow s})'})
%title({'P_{s\rightarrow s}'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})

% Inferred proportion female among primary
% res_arr(:,:,10:13) = [Q_pn Q_pp Q_sn Q_sp];
subplot (3,2,3)
temp = d2.res_arr(:,:,11)./(d2.res_arr(:,:,10) + d2.res_arr(:,:,11));
contour(d1.rp_arr,d1.p_obs_arr,temp','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Inferred proportion of','primary cases that are female'})
ylabel({'Observation probability (P_{obs})'},'FontSize',16)
%xlabel({'              Primary cases per cluster (R_p)'})

% Inferred proportion female among secondary
% res_arr(:,:,10:13) = [Q_pn Q_pp Q_sn Q_sp];
subplot (3,2,4)
temp = d2.res_arr(:,:,13)./(d2.res_arr(:,:,12) + d2.res_arr(:,:,13))
contour(d1.rp_arr,d1.p_obs_arr,temp','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Inferred proportion of','secondary cases that are female'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel({'              Primary cases per cluster (R_p)'})

% Inferred r_effective
subplot (3,2,5)
contour(d1.rp_arr,d1.p_obs_arr,R_inf,'ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Effective reproduction number'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel('Primaries per cluster (R_p)')

% Inferred odds ratio
% res_arr(:,:,10:13) = [Q_pn Q_pp Q_sn Q_sp];
subplot (3,2,6)
temp = d2.res_arr(:,:,10).*d2.res_arr(:,:,13)./(d2.res_arr(:,:,11) .* d2.res_arr(:,:,12))
% abs(temp) used because of some convergence difficulties yielding large negative
% values of OR for high p_obs and r_p
% contour(d1.rp_arr,d1.p_obs_arr,abs(temp)',[0.01 0.02 0.03 0.04 0.05 0.06],'ShowText','On','LineWidth',2)
contour(d1.rp_arr,d1.p_obs_arr,temp',[0.01 0.02 0.03 0.04 0.05 0.06],'ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Inferred odds ratio for','secondary cases being female'})
%ylabel({'Observation_ ','probability (P_{obs})'})
%xlabel('Primaries per cluster (R_p)')
text(.85,0.07,0,'Primary cases per cluster (R_p)','FontSize', 16)

% savepdf('Figures/040519_mpx_app')
exportgraphics(gcf, 'Figures/050125_nipah_app.pdf', 'Resolution', 300)

%%
% See script_disp050818 for the following figures:
%
% Figure for showing what the r_e is for various combinations of community
% and hospital transmission
%
% Goal: For a few R_eff, plot reduction in incidence as a function of
% decreased introduction and decreased transmission

%%
% See script_disp122018 for the following figures:
%
% Comparing MPX R in 1980s ot 2000s
%
% Analyzing MRSA data from Coll et al.
%% Updated figures for submission June 22, 2026

%% Figure: Bias in estimating the reproduction number

d1 = load('Data/011119_rp_pobs');
figure(1),clf
tiledlayout(2,1,'TileSpacing','tight','Padding','compact')
ax1 = nexttile;
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,9))','ShowText','On','LineWidth',2);
hold on
plot(min(d1.rp_arr), max(d1.p_obs_arr), 'r.', 'MarkerSize', 15)
hold off
set(gca, 'Layer', 'bottom')
axis square
set(gca,'FontSize',9)
title({'Bias of R_s estimate (\delta_R)'})
set(gca,'XTickLabel','')
tt = text(2.01, .9, {'R_s = 0.2'});
set(tt, 'FontWeight','bold')

ax2 = nexttile;
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,2,1,:,9))', 'ShowText','On','LineWidth',2);
hold on
plot(min(d1.rp_arr), max(d1.p_obs_arr), 'r.', 'MarkerSize', 15)
hold off
set(gca, 'Layer', 'bottom')
axis square
set(gca,'FontSize',9)
xlabel({'Primary cases per cluster (R_p)'})
tt = text(2.01, .9, {'R_s = 0.8'});
set(tt, 'FontWeight','bold')
text(0.82,.75,0,'Observation probability (p_{obs})','Rotation', 90)

clim(ax1, [-0.2 0.35])
clim(ax2, [-0.2 0.35])

set(gcf, 'position', [20, 20, 1400, 600])
exportgraphics(gcf, 'Figures/R_est_062226.pdf')

%% Figure: Attenuation of the observed odds ratio
figure(1),clf
d = load('Data/040419_ORx_example');

nn = 0;
for rr = 1:length(d.rs_arr)
    for ii = 1:length(d.inc_arr)
        nn = nn+1;
        subplot(length(d.rs_arr),length(d.inc_arr),nn)
        contour(d.rp_arr,d.p_obs_arr,squeeze(d.res_arr(:,:,rr,ii,2)'), [1.5 2 2.5 3 3.5],'ShowText','On','LineWidth',2)
        hold on;
        
        title(strcat(['R_s =',' ',num2str(d.rs_arr(rr)),'   I_p =',' ',num2str(d.inc_arr(ii,1))]));
        clim([1.5 3.5])
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
        hold on
        plot(1, 1, 'r.', 'MarkerSize', 15)
        hold off
        set(gca, 'Layer', 'bottom')
        text(1.01,0.95,'4')
    end
end

set(gcf, 'position', [20, 20, 1000, 800])
exportgraphics(gcf, 'Figures/ORx_062626.pdf')

%% Figure: Applying our methods to data from mpox in the Democratic 
% Republic of the Congo, 2013–2017

figure(1),clf

clear
% cannot find this data
d2 = load('Data/082225_2010s_mpox_OR.mat');
d1 = d2.d_mpox;
R_inf = squeeze(d1.rs_inf_res)';

% Prob observed primary is a primary
subplot (2,2,1)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,5)',[0.7, 0.8, 0.9, 1],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a case classified_ ','as primary is primary (P_{p\rightarrow p})'})

% Prob observed secondary is a secondary
subplot (2,2,2)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,6)','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a case classified_ ','as secondary is secondary (P_{s\rightarrow s})'})

% Inferred r_effective
subplot (2,2,3)
contour(d1.rp_arr,d1.p_obs_arr,R_inf,'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Effective reproduction number'})
xlabel(' ')
text(1.42,0.08,0,'Primary cases per cluster (R_p)','FontSize', 11)
ylabel(' ')
text(0.92, 0.85, 0, 'Observation probability (P_{obs})', 'FontSize', 11, 'Rotation', 90)

% Inferred odds ratio
subplot (2,2,4)
temp = d2.res_arr(:,:,11).*d2.res_arr(:,:,12)./(d2.res_arr(:,:,10) .* d2.res_arr(:,:,13));
contour(d1.rp_arr,d1.p_obs_arr,temp',[1.5, 1.55, 1.6, 1.65, 1.7],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Inferred odds ratio for primary','cases reporting animal contact'})

set(gcf, 'position', [20, 20, 800, 600])
exportgraphics(gcf, 'Figures/mpox_2010s_062626.pdf')

%% Figure: Probability that a true primary case is classified as a primary 
% or secondary case

d1 = load('Data/011119_rp_pobs');

figure(1),clf
set(gcf, 'Position', [10, 10, 1200, 900])
subplot(2,2,1)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,1))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that a primary case is_ ', 'classified as a primary case (C_{p\rightarrow p})'})
ylabel({'Observation_ ','probability (P_{obs})'},'FontSize',16)
subplot(2,2,2)
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,2))','ShowText','On','LineWidth',2)
set(gca,'FontSize',12)
title({'Probability that a primary case is_ ', 'classified as a secondary case (C_{p\rightarrow s})'})
set(gca,'YTickLabel','')
text(.5,-.05,0,'Primary cases per cluster (R_p)','FontSize', 16)

exportgraphics(gcf, 'Figures/class_prim_062626.pdf')

%% Figure: Probability that a true secondary case is classified as a 
% primary or secondary case

d1 = load('Data/011119_rp_pobs');
figure(1),clf
tiledlayout(3,2,'TileSpacing','loose','Padding','compact')

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,2,:,3))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a secondary case is_ ', 'classified as a primary case (C_{s\rightarrow p})'})
set(gca,'XTickLabel','')
tt = text(2.01, .9, {'Homogeneous','R_s = 0.2'});
set(tt,'FontWeight','bold')

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,2,:,4))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
set(gca,'XTickLabel','')
title({'Probability that a secondary case is_ ', 'classified as a secondary case (C_{s\rightarrow s})'})
set(gca,'YTickLabel','')


nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,3))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
ylabel({'Observation probability (P_{obs})'},'FontSize',12)
set(gca,'XTickLabel','')
tt = text(2.01, .9, {'Heterogeneous','R_s = 0.2'});
set(tt, 'FontWeight','bold')

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,4))','ShowText','On','LineWidth',2)
set(gca,'XTickLabel','')
set(gca,'FontSize',9)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')


nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,2,1,:,3))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
tt = text(2.01, .9, {'Heterogeneous','R_s = 0.8'});
set(tt, 'FontWeight','bold')

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,2,1,:,4))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
set(gca,'YTickLabel','')
text(.56,-.08,0,'Primary cases per cluster (R_p)', 'FontSize',12)

set(gcf, 'position', [20, 20, 1000, 900])
exportgraphics(gcf, 'Figures/class_sec_062226.pdf')


%% Figure: Probabilities that observed cases are correctly classified as 
% primary or secondary cases

d1 = load('Data/011119_rp_pobs');
figure(1),clf
tiledlayout(3,4,'TileSpacing','loose','Padding','compact')

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,2,:,5))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a case_ ','classified as a primary case_ ', 'is a primary case (P_{p\rightarrow p})'})
set(gca,'XTickLabel','')

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,2,:,6))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
set(gca,'XTickLabel','')
title({'Probability that a case_ ','classified as a secondary case_ ', 'is a secondary case (P_{s\rightarrow s})'})
set(gca,'YTickLabel','')

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,2,:,7))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
title({'Probability that an observed_ ','case is classified correctly (\Theta)'})

nexttile
tt = text(-0.25, 0.9, {'Homogeneous','R_s = 0.2'});
set(tt, 'FontSize', 10, 'FontWeight','bold')
axis off


nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,5))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
set(gca,'XTickLabel','')
ylabel({'Observation probability (P_{obs})'})

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,6))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,1,1,:,7))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')

nexttile
tt = text(-0.25, 0.9, {'Heterogeneous','R_s = 0.2'});
set(tt, 'FontSize', 10, 'FontWeight','bold')
axis off


nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,2,1,:,5))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,2,1,:,6))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
xlabel({'Primary cases per cluster (R_p)'})
set(gca,'YTickLabel','')

nexttile
contour(d1.rp_arr,d1.p_obs_arr,squeeze(d1.res_arr(:,2,1,:,7))','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
set(gca,'YTickLabel','')

nexttile
tt = text(-0.25, 0.9, {'Heterogeneous','R_s = 0.8'});
set(tt, 'FontSize', 10, 'FontWeight','bold')
axis off

set(gcf, 'position', [20, 20, 1200, 1000])
exportgraphics(gcf, 'Figures/class_accuracy_062226.pdf')


%% Figure: Mpox in the Democratic Republic of the Congo, 1981-86: Accuracy 
% of classifying cases as primary versus secondary, and implications for 
% the gender distribution of cases

figure(1),clf

clear
d2 = load('Data/111925_1980s_mpox_OR.mat');
d1 = d2.d_mpx;
R_inf = squeeze(d1.rs_inf_res)';

% Prob observed primary is a primary
subplot (2,2,1)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,5)',[0.75, 0.8, 0.85, 0.9, 0.95, 1],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a case classified_ ','as primary is primary (P_{p\rightarrow p})'})

% Prob observed secondary is a secondary
subplot (2,2,2)
contour(d1.rp_arr,d1.p_obs_arr,d2.res_arr(:,:,6)','ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Probability that a case classified_ ','as secondary is secondary (P_{s\rightarrow s})'})

% Inferred r_effective
subplot (2,2,3)
contour(d1.rp_arr,d1.p_obs_arr,R_inf,'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Effective reproduction number'})
xlabel(' ')
text(1.25,0.08,0,'Primary cases per cluster (R_p)','FontSize', 11)
ylabel(' ')
text(0.94, 0.85, 0, 'Observation probability (P_{obs})', 'FontSize', 11, 'Rotation', 90)

% Inferred odds ratio
subplot (2,2,4)
temp = d2.res_arr(:,:,13).*d2.res_arr(:,:,10)./(d2.res_arr(:,:,12) .* d2.res_arr(:,:,11));
contour(d1.rp_arr,d1.p_obs_arr,abs(temp)',[1.9, 2.1, 2.3, 2.5, 3, 5],'ShowText','On','LineWidth',2)
set(gca,'FontSize',9)
title({'Inferred odds ratio for secondary','cases being female'})

set(gcf, 'position', [20, 20, 800, 600])
exportgraphics(gcf, 'Figures/mpox_1980s_062226.pdf')


%% Eval 3: Role of Prior on MIDAS Component --- Full Sample %%

clear;
clear all

%% Load Data
%cd '/Users/dk/triton_work/MIDAS/Bayesian-MIDAS'

addpath('../Matlab')
% Load Data
load('../Data/UK_dat_2024.mat');

rng(1);

%% GIGG Models
% Trend-SV-GIGG
mod1 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_trend_sv_almon_ortho_groupsparse.mat');
mod1 = mod1.output;
mod1.resid_all([1:5 12 end],:) = [];
mod1.crps_all([1:5 12 end],:) = [];

% Trend-SV-GIGG No Sparsification
mod2 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_trend_sv_almon_ortho.mat');
mod2 = mod2.output;
mod2.resid_all([1:5 12 end],:) = [];
mod2.crps_all([1:5 12 end],:) = [];

% Trend-SV-HS
mod3 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_trend_sv_umidas_ortho_groupsparse.mat');
mod3 = mod3.output;
mod3.resid_all([1:5 12 end],:) = [];
mod3.crps_all([1:5 12 end],:) = [];

% Trend-SV-SS
mod4 = load('Output_iterated/results_iteratednowcasts_hs_oldsv_correcttrend_newdat_trend_sv_almon_groupsparse.mat');
mod4 = mod4.output;
mod4.resid_all([1:5 12 end],:) = [];
mod4.crps_all([1:5 12 end],:) = [];

% UMIDAS Combination
mod5 = load('Output_iterated/results_iteratednowcasts_ssvs_oldsv_correcttrend_newdat_trend_sv_almon_groupsparse.mat');
mod5 = mod5.output;
mod5.resid_all([1:5 12 end],:) = [];
mod5.crps_all([1:5 12 end],:) = [];

% Trend-SV-GIGG UMIDAS
mod6 = load('Output_iterated/results_iteratednowcasts_midascomb_umidas.mat');
mod6 = mod6.output;
mod6.resid_all([1:5 12 end],:) = [];
mod6.crps_all([1:5 12 end],:) = [];


yf = yf_all;

res_mod1 = mod1.resid_all;
res_mod2 = mod2.resid_all;
res_mod3 = mod3.resid_all;
res_mod4 = mod4.resid_all;
res_mod5 = mod5.resid_all;
res_mod6 = mod6.resid_all;

%% Metrics including all data points

    % GFC
rmsfe_gfc_mod1 = std(res_mod1(:,1:12)')';
rmsfe_gfc_mod2 = std(res_mod2(:,1:12)')';
rmsfe_gfc_mod3 = std(res_mod3(:,1:12)')';
rmsfe_gfc_mod4 = std(res_mod4(:,1:12)')';
rmsfe_gfc_mod5 = std(res_mod5(:,1:12)')';
rmsfe_gfc_mod6 = std(res_mod6(:,1:12)')';


rtcrps_gfc_mod1= mean(mod1.crps_all(:,1:12),2);
rtcrps_gfc_mod2= mean(mod2.crps_all(:,1:12),2);
rtcrps_gfc_mod3= mean(mod3.crps_all(:,1:12),2);
rtcrps_gfc_mod4= mean(mod4.crps_all(:,1:12),2);
rtcrps_gfc_mod5= mean(mod5.crps_all(:,1:12),2);
rtcrps_gfc_mod6= mean(mod6.crps_all(:,1:12),2);


    % Tranquil
rmsfe_tranq_mod1 = std(res_mod1(:,13:52)')';
rmsfe_tranq_mod2 = std(res_mod2(:,13:52)')';
rmsfe_tranq_mod3 = std(res_mod3(:,13:52)')';
rmsfe_tranq_mod4 = std(res_mod4(:,13:52)')';
rmsfe_tranq_mod5 = std(res_mod5(:,13:52)')';
rmsfe_tranq_mod6 = std(res_mod6(:,13:52)')';


rtcrps_tranq_mod1= mean(mod1.crps_all(:,13:52),2);
rtcrps_tranq_mod2= mean(mod2.crps_all(:,13:52),2);
rtcrps_tranq_mod3= mean(mod3.crps_all(:,13:52),2);
rtcrps_tranq_mod4= mean(mod4.crps_all(:,13:52),2);
rtcrps_tranq_mod5= mean(mod5.crps_all(:,13:52),2);
rtcrps_tranq_mod6= mean(mod6.crps_all(:,13:52),2);


    % Pandemic
rmsfe_pandemic_mod1 = std(res_mod1(:,53:end)')';
rmsfe_pandemic_mod2 = std(res_mod2(:,53:end)')';
rmsfe_pandemic_mod3 = std(res_mod3(:,53:end)')';
rmsfe_pandemic_mod4 = std(res_mod4(:,53:end)')';
rmsfe_pandemic_mod5 = std(res_mod5(:,53:end)')';
rmsfe_pandemic_mod6 = std(res_mod6(:,53:end)')';


rtcrps_pandemic_mod1= mean(mod1.crps_all(:,53:end),2);
rtcrps_pandemic_mod2= mean(mod2.crps_all(:,53:end),2);
rtcrps_pandemic_mod3= mean(mod3.crps_all(:,53:end),2);
rtcrps_pandemic_mod4= mean(mod4.crps_all(:,53:end),2);
rtcrps_pandemic_mod5= mean(mod5.crps_all(:,53:end),2);
rtcrps_pandemic_mod6= mean(mod6.crps_all(:,53:end),2);



    % Including the Pandemic
rmsfe_post_mod1 = std(res_mod1(:,1:end)')';
rmsfe_post_mod2 = std(res_mod2(:,1:end)')';
rmsfe_post_mod3 = std(res_mod3(:,1:end)')';
rmsfe_post_mod4 = std(res_mod4(:,1:end)')';
rmsfe_post_mod5 = std(res_mod5(:,1:end)')';
rmsfe_post_mod6 = std(res_mod6(:,1:end)')';


rtcrps_post_mod1= mean(mod1.crps_all(:,1:end),2);
rtcrps_post_mod2= mean(mod2.crps_all(:,1:end),2);
rtcrps_post_mod3= mean(mod3.crps_all(:,1:end),2);
rtcrps_post_mod4= mean(mod4.crps_all(:,1:end),2);
rtcrps_post_mod5= mean(mod5.crps_all(:,1:end),2);
rtcrps_post_mod6= mean(mod6.crps_all(:,1:end),2);


%% Figure 1: Eval Graph
    % Upper two panels are RMSFE (pre and with pandemic)
    % Lower two panels are CRPS (pre and with pandemic)

fig1 = figure;

 subplot(2,4,1);
plot(rmsfe_post_mod1(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="-",MarkerSize=2)
hold on
plot(rmsfe_post_mod2(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="--",MarkerSize=2)
plot(rmsfe_post_mod6(1:end),'LineWidth',2,Color="#A73E73",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_post_mod3(1:end),'LineWidth',2,Color="#3EA772",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_post_mod4(1:end),'LineWidth',2,Color="#8F968B",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_post_mod5(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19]) 
xlim([0 20])
ylim([0 4.5])
h=get(fig1,'CurrentAxes')
ylabel('RMSFE','FontSize',16)
set(h, 'FontSize', 12) 
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
title('Full Sample','FontSize',16)

subplot(2,4,2);
plot(rmsfe_gfc_mod1(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="-",MarkerSize=2)
hold on
plot(rmsfe_gfc_mod2(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="--",MarkerSize=2)
plot(rmsfe_gfc_mod6(1:end),'LineWidth',2,Color="#A73E73",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_gfc_mod3(1:end),'LineWidth',2,Color="#3EA772",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_gfc_mod4(1:end),'LineWidth',2,Color="#8F968B",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_gfc_mod5(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19]) 
xlim([0 20])
ylim([0 4.5])
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
title('GFC','FontSize',16)

subplot(2,4,3);
plot(rmsfe_tranq_mod1(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="-",MarkerSize=2)
hold on
plot(rmsfe_tranq_mod2(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="--",MarkerSize=2)
plot(rmsfe_tranq_mod6(1:end),'LineWidth',2,Color="#A73E73",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_tranq_mod3(1:end),'LineWidth',2,Color="#3EA772",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_tranq_mod4(1:end),'LineWidth',2,Color="#8F968B",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_tranq_mod5(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19]) 
xlim([0 20])
ylim([0 4.5])
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
title('Tranquil','FontSize',16)

subplot(2,4,4);
plot(rmsfe_pandemic_mod1(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="-",MarkerSize=2)
hold on
plot(rmsfe_pandemic_mod2(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="--",MarkerSize=2)
plot(rmsfe_pandemic_mod6(1:end),'LineWidth',2,Color="#A73E73",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_pandemic_mod3(1:end),'LineWidth',2,Color="#3EA772",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_pandemic_mod4(1:end),'LineWidth',2,Color="#8F968B",Marker="o",LineStyle="-",MarkerSize=2)
plot(rmsfe_pandemic_mod5(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19]) 
xlim([0 20])
ylim([0 10])
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
title('Pandemic','FontSize',16)

subplot(2,4,5);
plot(rtcrps_post_mod1(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="-",MarkerSize=2)
hold on
plot(rtcrps_post_mod2(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="--",MarkerSize=2)
plot(rtcrps_post_mod6(1:end),'LineWidth',2,Color="#A73E73",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_post_mod3(1:end),'LineWidth',2,Color="#3EA772",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_post_mod4(1:end),'LineWidth',2,Color="#8F968B",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_post_mod5(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19]) 
xlim([0 20])
ylim([0 2])
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 
ylabel('CRPS','FontSize',16)

subplot(2,4,6);
plot(rtcrps_gfc_mod1(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="-",MarkerSize=2)
hold on
plot(rtcrps_gfc_mod2(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="--",MarkerSize=2)
plot(rtcrps_gfc_mod6(1:end),'LineWidth',2,Color="#A73E73",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_gfc_mod3(1:end),'LineWidth',2,Color="#3EA772",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_gfc_mod4(1:end),'LineWidth',2,Color="#8F968B",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_gfc_mod5(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19]) 
xlim([0 20])
ylim([0 2])
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 

subplot(2,4,7);
plot(rtcrps_tranq_mod1(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="-",MarkerSize=2)
hold on
plot(rtcrps_tranq_mod2(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="--",MarkerSize=2)
plot(rtcrps_tranq_mod6(1:end),'LineWidth',2,Color="#A73E73",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_tranq_mod3(1:end),'LineWidth',2,Color="#3EA772",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_tranq_mod4(1:end),'LineWidth',2,Color="#8F968B",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_tranq_mod5(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19]) 
xlim([0 20])
ylim([0 2])
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 

subplot(2,4,8);
plot(rtcrps_pandemic_mod1(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="-",MarkerSize=2)
hold on
plot(rtcrps_pandemic_mod2(1:end),'LineWidth',2,Color="#000000",Marker="o",LineStyle="--",MarkerSize=2)
plot(rtcrps_pandemic_mod6(1:end),'LineWidth',2,Color="#A73E73",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_pandemic_mod3(1:end),'LineWidth',2,Color="#3EA772",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_pandemic_mod4(1:end),'LineWidth',2,Color="#8F968B",Marker="o",LineStyle="-",MarkerSize=2)
plot(rtcrps_pandemic_mod5(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19]) 
xlim([0 20])
ylim([0 6])
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 

hL =  legend('GIGG','GIGG No Spars.','U-MIDAS-Comb.','GIGG U-MIDAS','HS','SSVS','NumColumns',3);
newPosition = [0.6 0.513 0.01 0.01] ;
set(hL,'Position', newPosition,'FontSize',14);
h=fig1;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
%set(h,'PaperPosition', [0 0 1 1]);
set(h,'PaperPosition', [-0.05 0 1.1 1]);

% Give common xlabel, ylabel and title to your figure
han=axes(h,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Days Until GDP Release','FontSize',16);

xh = get(han,'xlabel') % handle to the label object
p = get(xh,'position') % get the current position property
p(2) = 1.3*p(2)         % double the distance, 
set(xh,'position',p)   % set the new position


figname = ['Figures/eval_fig4_fullsamp_iterated.pdf'];
saveas(han,figname)
saveas(han,['Figures/eval_fig4_fullsamp_iterated.fig'])




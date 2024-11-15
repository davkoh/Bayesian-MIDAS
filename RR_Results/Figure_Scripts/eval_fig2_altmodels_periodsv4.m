%% Eval 2: Comparison against alternative Models --- Full Sample %%

clear;
clear all

%% Load Data
%cd '/Users/dk/triton_work/MIDAS/Bayesian-MIDAS'

addpath('../Matlab')
% Load Data
load('../Data/UK_dat_2024.mat');

rng(1);

% Trend-SV-GIGG
mod1 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_trend_sv_almon_ortho_groupsparse.mat');
mod1 = mod1.output;
mod1.resid_all([1:5 12 end],:) = [];
mod1.crps_all([1:5 12 end],:) = [];

% SV-GIGG 
mod2 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_sv_almon_ortho_groupsparse.mat');
mod2 = mod2.output;
mod2.resid_all([1:5 12 end],:) = [];
mod2.crps_all([1:5 12 end],:) = [];

% GIGG
mod3 = load('Output_iterated/results_iteratednowcasts_oldgigg_plain2_bg_0.5_almon_ortho_groupsparse.mat');
mod3 = mod3.output;
mod3.resid_all([1:5 12 end],:) = [];
mod3.crps_all([1:5 12 end],:) = [];

% HS
mod4 = load('Output_iterated/results_iteratednowcasts_hs_oldsv_correcttrend_newdat_almon_groupsparse.mat');
mod4 = mod4.output;
mod4.resid_all([1:5 12 end],:) = [];
mod4.crps_all([1:5 12 end],:) = [];

% M&S
mod5 = load('Output_iterated/results_iteratednowcasts_ms2_newdat.mat');
mod5 = mod5.output;
mod5.resid_all([1:5 12 end],:) = [];
mod5.crps_all([1:5 12 end],:) = [];

% PCA-MIDAS
mod6 = load('Output_iterated/results_iteratednowcasts_pca_newsv_newdat_umidas.mat');
mod6 = mod6.output;
mod6.resid_all([1:5 12 end],:) = [];
mod6.crps_all([1:5 12 end],:) = [];

% Combined MIDAS
mod7 = load('Output_iterated/results_iteratednowcasts_midascomb_umidas.mat');
mod7 = mod7.output;
mod7.resid_all([1:5 12 end],:) = [];
mod7.crps_all([1:5 12 end],:) = [];



yf = yf_all;

res_mod1 = mod1.resid_all;
res_mod2 = mod2.resid_all;
res_mod3 = mod3.resid_all;
res_mod4 = mod4.resid_all;
res_mod5 = mod5.resid_all;
res_mod6 = mod6.resid_all;
res_mod7 = mod7.resid_all;

%% Metrics GFC

rmsfe_gfc_mod1 = std(res_mod1(:,1:12)')';
rmsfe_gfc_mod2 = std(res_mod2(:,1:12)')';
rmsfe_gfc_mod3 = std(res_mod3(:,1:12)')';
rmsfe_gfc_mod4 = std(res_mod4(:,1:12)')';
rmsfe_gfc_mod5 = std(res_mod5(:,1:12)')';
rmsfe_gfc_mod6 = std(res_mod6(:,1:12)')';
rmsfe_gfc_mod7 = std(res_mod7(:,1:12)')';


rtcrps_gfc_mod1= mean(mod1.crps_all(:,1:12),2);
rtcrps_gfc_mod2= mean(mod2.crps_all(:,1:12),2);
rtcrps_gfc_mod3= mean(mod3.crps_all(:,1:12),2);
rtcrps_gfc_mod4= mean(mod4.crps_all(:,1:12),2);
rtcrps_gfc_mod5= mean(mod5.crps_all(:,1:12),2);
rtcrps_gfc_mod6= mean(mod6.crps_all(:,1:12),2);
rtcrps_gfc_mod7= mean(mod7.crps_all(:,1:12),2);

%% Metrics Tranquil

rmsfe_tranq_mod1 = std(res_mod1(:,13:52)')';
rmsfe_tranq_mod2 = std(res_mod2(:,13:52)')';
rmsfe_tranq_mod3 = std(res_mod3(:,13:52)')';
rmsfe_tranq_mod4 = std(res_mod4(:,13:52)')';
rmsfe_tranq_mod5 = std(res_mod5(:,13:52)')';
rmsfe_tranq_mod6 = std(res_mod6(:,13:52)')';
rmsfe_tranq_mod7 = std(res_mod7(:,13:52)')';


rtcrps_tranq_mod1= mean(mod1.crps_all(:,13:52),2);
rtcrps_tranq_mod2= mean(mod2.crps_all(:,13:52),2);
rtcrps_tranq_mod3= mean(mod3.crps_all(:,13:52),2);
rtcrps_tranq_mod5= mean(mod5.crps_all(:,13:52),2);
rtcrps_tranq_mod6= mean(mod6.crps_all(:,13:52),2);
rtcrps_tranq_mod4= mean(mod4.crps_all(:,13:52),2);
rtcrps_tranq_mod7= mean(mod7.crps_all(:,13:52),2);

%% Metrics Covid

rmsfe_covid_mod1 = std(res_mod1(:,53:end)')';
rmsfe_covid_mod2 = std(res_mod2(:,53:end)')';
rmsfe_covid_mod3 = std(res_mod3(:,53:end)')';
rmsfe_covid_mod4 = std(res_mod4(:,53:end)')';
rmsfe_covid_mod5 = std(res_mod5(:,53:end)')';
rmsfe_covid_mod6 = std(res_mod6(:,53:end)')';
rmsfe_covid_mod7 = std(res_mod7(:,53:end)')';


rtcrps_covid_mod1= mean(mod1.crps_all(:,53:end),2);
rtcrps_covid_mod2= mean(mod2.crps_all(:,53:end),2);
rtcrps_covid_mod3= mean(mod3.crps_all(:,53:end),2);
rtcrps_covid_mod4= mean(mod4.crps_all(:,53:end),2);
rtcrps_covid_mod5= mean(mod5.crps_all(:,53:end),2);
rtcrps_covid_mod6= mean(mod6.crps_all(:,53:end),2);
rtcrps_covid_mod7= mean(mod7.crps_all(:,53:end),2);


%% Figure 1: Eval Graph

fig1 = figure;


subplot(2,3,1);
plot(rmsfe_gfc_mod1(1:end),'LineWidth',3,Color="#000000",LineStyle="-")
hold on
plot(rmsfe_gfc_mod2(1:end),'LineWidth',2,Color="#000000",LineStyle="--")
plot(rmsfe_gfc_mod3(1:end),'LineWidth',2,Color="#000000",LineStyle=":")
plot(rmsfe_gfc_mod4(1:end),'LineWidth',2,Color = "#8F968B",Marker="square",LineStyle="-",MarkerSize=2)
plot(rmsfe_gfc_mod5(1:end),'LineWidth',2,Color="#3EA772",Marker="+",LineStyle="-",MarkerSize=2)
plot(rmsfe_gfc_mod6(1:end),'LineWidth',2,Color="#ffa79b",Marker="square",LineStyle="--",MarkerSize=2)
plot(rmsfe_gfc_mod7(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19])
xlim([0 20])
ylim([0 1.5])
h=get(fig1,'CurrentAxes')
ylabel('RMSFE','FontSize',16)
set(h, 'FontSize', 12) 
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
title('GFC','FontSize',16)

 subplot(2,3,2);
plot(rmsfe_tranq_mod1(1:end),'LineWidth',3,Color="#000000",LineStyle="-")
hold on
plot(rmsfe_tranq_mod2(1:end),'LineWidth',2,Color="#000000",LineStyle="--")
plot(rmsfe_tranq_mod3(1:end),'LineWidth',2,Color="#000000",LineStyle=":")
plot(rmsfe_tranq_mod4(1:end),'LineWidth',2,Color = "#8F968B",Marker="square",LineStyle="-",MarkerSize=2)
plot(rmsfe_tranq_mod5(1:end),'LineWidth',2,Color="#3EA772",Marker="+",LineStyle="-",MarkerSize=2)
plot(rmsfe_tranq_mod6(1:end),'LineWidth',2,Color="#ffa79b",Marker="square",LineStyle="--",MarkerSize=2)
plot(rmsfe_tranq_mod7(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19])
xlim([0 20])
%ylim([0 1.1])
ylim([0 1.5])
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
title('Tranquil','FontSize',16)

subplot(2,3,3);
plot(rmsfe_covid_mod1(1:end),'LineWidth',3,Color="#000000",LineStyle="-")
hold on
plot(rmsfe_covid_mod2(1:end),'LineWidth',2,Color="#000000",LineStyle="--")
plot(rmsfe_covid_mod3(1:end),'LineWidth',2,Color="#000000",LineStyle=":")
plot(rmsfe_covid_mod4(1:end),'LineWidth',2,Color = "#8F968B",Marker="square",LineStyle="-",MarkerSize=2)
plot(rmsfe_covid_mod5(1:end),'LineWidth',2,Color="#3EA772",Marker="+",LineStyle="-",MarkerSize=2)
plot(rmsfe_covid_mod6(1:end),'LineWidth',2,Color="#ffa79b",Marker="square",LineStyle="--",MarkerSize=2)
plot(rmsfe_covid_mod7(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19])
%ylim([0 6])
xlim([0 20])
%ylim([0 1.1])
ylim([0 10])
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 
title('Pandemic','FontSize',16)


subplot(2,3,4);
plot(rtcrps_gfc_mod1(1:end),'LineWidth',3,Color="#000000",LineStyle="-")
hold on
plot(rtcrps_gfc_mod2(1:end),'LineWidth',2,Color="#000000",LineStyle="--")
plot(rtcrps_gfc_mod3(1:end),'LineWidth',2,Color="#000000",LineStyle=":")
plot(rtcrps_gfc_mod4(1:end),'LineWidth',2,Color = "#8F968B",Marker="square",LineStyle="-",MarkerSize=2)
plot(rtcrps_gfc_mod5(1:end),'LineWidth',2,Color="#3EA772",Marker="+",LineStyle="-",MarkerSize=2)
plot(rtcrps_gfc_mod6(1:end),'LineWidth',2,Color="#ffa79b",Marker="square",LineStyle="--",MarkerSize=2)
plot(rtcrps_gfc_mod7(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19])
xlim([0 20])
ylim([0 1.5])
h=get(fig1,'CurrentAxes')
ylabel('CRPS','FontSize',16)
set(h, 'FontSize', 12) 
xticklabels({'135','120','110','95','85','75','60','50','35','15'})

subplot(2,3,5);
plot(rtcrps_tranq_mod1(1:end),'LineWidth',3,Color="#000000",LineStyle="-")
hold on
plot(rtcrps_tranq_mod2(1:end),'LineWidth',2,Color="#000000",LineStyle="--")
plot(rtcrps_tranq_mod3(1:end),'LineWidth',2,Color="#000000",LineStyle=":")
plot(rtcrps_tranq_mod4(1:end),'LineWidth',2,Color = "#8F968B",Marker="square",LineStyle="-",MarkerSize=2)
plot(rtcrps_tranq_mod5(1:end),'LineWidth',2,Color="#3EA772",Marker="+",LineStyle="-",MarkerSize=2)
plot(rtcrps_tranq_mod6(1:end),'LineWidth',2,Color="#ffa79b",Marker="square",LineStyle="--",MarkerSize=2)
plot(rtcrps_tranq_mod7(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19])
%ylim([0 6])
xlim([0 20])
ylim([0 1.5])
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 

subplot(2,3,6);
plot(rtcrps_covid_mod1(1:end),'LineWidth',3,Color="#000000",LineStyle="-")
hold on
plot(rtcrps_covid_mod2(1:end),'LineWidth',2,Color="#000000",LineStyle="--")
plot(rtcrps_covid_mod3(1:end),'LineWidth',2,Color="#000000",LineStyle=":")
plot(rtcrps_covid_mod4(1:end),'LineWidth',2,Color = "#8F968B",Marker="square",LineStyle="-",MarkerSize=2)
plot(rtcrps_covid_mod5(1:end),'LineWidth',2,Color="#3EA772",Marker="+",LineStyle="-",MarkerSize=2)
plot(rtcrps_covid_mod6(1:end),'LineWidth',2,Color="#ffa79b",Marker="square",LineStyle="--",MarkerSize=2)
plot(rtcrps_covid_mod7(1:end),'LineWidth',2,Color="#FFD99B",Marker="o",LineStyle="-",MarkerSize=2)
xticks([1,3,5,7,9,11,13,15,17,19])
ylim([0 7])
xlim([0 20])
xticklabels({'135','120','110','95','85','75','60','50','35','15'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 12) 

hL =  legend('T-SV-BMIDAS (GIGG)','SV-BMIDAS (GIGG)' ,'BMIDAS (GIGG)','BMIDAS (HS)','MS2021','U-BMIDAS-PCA','U-BMIDAS-Comb','NumColumns',4);
newPosition = [0.5 0.49 0.01 0.01] ;
set(hL,'Position', newPosition,'FontSize',14);
h=fig1;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [-0.05 0 1.1 1]);
%set(h,'PaperPosition', [0 0 1 1]);

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

figname = ['Figures/eval_fig2_periods_iterated.pdf'];
saveas(han,figname)




%% Eval Graphs MF treatment + UMIDAS %%
% V2: now also includes UMIDAS + new colours and point markers


clear;
clear all

%% Load Data
cd 'C:\Users\dkohn\OneDrive - Heriot-Watt University\Education\BoE Nowcasting Project 2022\Estimation'

addpath('Data')
addpath('Mfiles')
addpath('D:\Downloads\GIGG-Nowcast-Output')
addpath('C:\Users\dkohn\OneDrive - Heriot-Watt University\Education\BoE Nowcasting Project 2022\Estimation\Output')
addpath('C:\Users\dkohn\OneDrive - Heriot-Watt University\Education\BoE Nowcasting Project 2022\Functions')
%addpath('N:\CECD\10. Personal\David Kohns\Code\Application\Macro Application\Results\Updated Sample')
% Load Data
load('UKdata_11_11_1980_2021_GDPCONINVHOUR_sur_act_lab_mort_ie_vis.mat');


%% OG GIGG Models
% Trend-SV-t-Almon
mod1 = load('out_now_gigg_sv_t_redocheck4_0.0217390.5_almon_ortho_gsavs_final.mat');
mod1 = mod1.output;
mod1.rmsfe_all([1:6],:) = [];
mod1.crps_all([1:6],:) = [];
mod1.rmsfe_all([2 8 15 21 25],:) = [];
mod1.crps_all([[2 8 15 21 25]],:) = [];

%% GIGG Model HS priors on SW components
% No t-errors
mod2 = load("output_hst_sw_corrected_0.020833_0.5_almon_ortho_groupsparse.mat")
mod2 = mod2.output;
mod2.resid_all([1:5 end],:) = [];
mod2.crps_all([1:5 end],:) = [];

% No t-errors
mod3 = load("output_hs_sw_corrected_0.020833_0.5_almon_ortho_groupsparse.mat")
mod3 = mod3.output;
mod3.resid_all([1:5 end],:) = [];
mod3.crps_all([1:5 end],:) = [];

%% Retrieve Scores for Non-GIGG Models
nfor = 43; %%%%%%%%%%% Needs to change
yf = y_q(end-nfor+1:end,1);
yf_dfm = y_q(end-43+1:end,1);

res_mod1 = mod1.rmsfe_all;
res_mod2 = mod2.resid_all;
res_mod3 = mod3.resid_all;


%% GIGG Models
    % Pre Pandemic
rmsfe_pre_mod1 = std(res_mod1(:,1:36)')';
rmsfe_pre_mod2 = std(res_mod2(:,1:35)')';
rmsfe_pre_mod3 = std(res_mod3(:,1:35)')';

rtcrps_pre_mod1 = mean(mod1.crps_all(:,33:73),2);
rtcrps_pre_mod2 = mean(mod2.crps_all(:,1:35),2);
rtcrps_pre_mod3 = mean(mod3.crps_all(:,1:35),2);


    % Including the Pandemic
rmsfe_post_mod1 = std(res_mod1(:,1:end)')';
rmsfe_post_mod2 = std(res_mod2(:,1:end)')';
rmsfe_post_mod3 = std(res_mod3(:,1:end)')';

rtcrps_post_mod1= mean(mod1.crps_all(:,44:end),2);
rtcrps_post_mod2= mean(mod2.crps_all(:,1:end),2);
rtcrps_post_mod3= mean(mod3.crps_all(:,1:end),2);




%% Figure 1: Eval Graph
    % Upper two panels are RMSFE (pre and with pandemic)
    % Lower two panels are CRPS (pre and with pandemic)

fig1 = figure;

 subplot(2,2,1);
plot(rmsfe_pre_mod1(1:end),'LineWidth',1.5)
hold on
plot(rmsfe_pre_mod2(1:end),'LineWidth',1.5)
plot(rmsfe_pre_mod3(1:end),'LineWidth',1.5)
xticks([1,4,8,12,15,18])
ylim([0.15 0.45])
xlim([0 22])
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 15) 
xticklabels({'135','115','95','75','55','35'})
title('RMSFE Pre-Pandemic','FontSize',16)

subplot(2,2,2);
plot(rmsfe_post_mod1(1:end),'LineWidth',1.5)
hold on
plot(rmsfe_post_mod2(1:end),'LineWidth',1.5)
plot(rmsfe_post_mod3(1:end),'LineWidth',1.5)
xticks([1,4,8,12,15,18])
ylim([0.5 7])
xlim([0 22])
xticklabels({'135','115','95','75','55','35'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 15) 
title('RMSFE Including Pandemic','FontSize',16)

subplot(2,2,3);
plot(rtcrps_pre_mod1(1:end),'LineWidth',1.5)
hold on
plot(rtcrps_pre_mod2(1:end),'LineWidth',1.5)
plot(rtcrps_pre_mod3(1:end),'LineWidth',1.5)
xticks([1,4,8,12,15,18])
ylim([0.10 0.25])
xlim([0 22])
xticklabels({'135','115','95','75','55','35'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 15) 
title('CRPS Pre-Pandemic','FontSize',16)

subplot(2,2,4);
plot(rtcrps_post_mod1(1:end),'LineWidth',1.5)
hold on
plot(rtcrps_post_mod2(1:end),'LineWidth',1.5)
plot(rtcrps_post_mod3(1:end),'LineWidth',1.5)
xticks([1,4,8,12,15,18])
ylim([0 3])
xlim([0 22])
xticklabels({'135','115','95','75','55','35'})
h=get(fig1,'CurrentAxes')
set(h, 'FontSize', 15) 
title('CRPS Including Pandemic','FontSize',16)


hL =  legend('KP (2023)','Trend(HS)-SV(HS)-GIGG','Trend(HS)-SV(HS)-t-GIGG');
%set(fig1,'Position',get(0,'Screensize'))

% Deal with Legend & Readability
%newPosition = [0.05,0.53,0.15,0.098406745154051];
newPosition = [0.8 0.8 0.1 0.1] ;
set(hL,'Position', newPosition,'FontSize',14);
h=fig1;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

% Give common xlabel, ylabel and title to your figure
han=axes(h,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Days Until GDP Release','FontSize',16);
figname = ['D:\Github\Bayesian-MIDAS\Output/eval_nowcasting.pdf'];
saveas(han,figname)
saveas(han,['D:\Github\Bayesian-MIDAS\Output/eval_nowcasting.fig'])



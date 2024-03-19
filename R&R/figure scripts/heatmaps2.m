%% Inclusion probability maps with the new models

cd '/Users/dk/Library/CloudStorage/OneDrive-Heriot-WattUniversity/Education/BoE Nowcasting Project 2022/Scripts from PC/Downloads'

clear all

%% To Do's
% Inlcude in the a vector of the length of groups which indicates
% overarching macro groups
%   Could be done by saving the 'varnames_qm' vector and comparing that with
%   the groups in the real time calendar
% Include capabilities to handle lags

%% Preload an output structure whith quarters available
old = load('/Users/dk/Library/CloudStorage/OneDrive-Heriot-WattUniversity/Education/BoE Nowcasting Project 2022/Estimation/Output/output_bsts_sv_t_low_wlags_noinfl')

%% Load Actual Output structure
load('out_now_gigg_sv_t_redocheck4_0.0217390.5_almon_ortho_gsavs_final.mat')
almon = output;
almon.incl([1:6],:,:) = [];
almon.incl([2 8 15 21 25],:,:) = [];

load('output_hs_0.020833_0.5_almon_ortho_groupsparse.mat')
hs_sw_svt = output;
hs_sw_svt.incl([1:5 end],:,:) = [];

load('output_hier_0.020833_0.5_almon_ortho_groupsparse.mat')
hs_sw_sv = output;
hs_sw_sv.incl([1:5 end],:,:) = [];

load('varnames_qm.mat')

%% Creating graphs
nfor = length(output.yf);
yf = output.yf;
d_q = old.output.d_q(end-nfor+1:end);
startdate = datenum(d_q(1));
enddate = datenum(d_q(end));
forecastdate = d_q; 
y = output.y;
savedir = '/Users/dk/Documents/GitHub/Bayesian-MIDAS/Output';
modelname1 = 'heat_combined.pdf' % Combined Heatmap plots!!!
modelname2 = 'average_cumulative_prob_prepandemic.pdf' % Average Cum Prob Plots before the pandemic !!!
modelname3 = 'average_cumulative_prob_postpandemic.pdf' % Average Cum Prob Plots before the pandemic !!!

%load('UKdata_11_11_1980_2021_GDPCONINVHOUR_sur_act_lab_mort_ie_vis.mat');

%% Find Grouping
T2 = readtable('/Users/dk/Library/CloudStorage/OneDrive-Heriot-WattUniversity/Education/BoE Nowcasting Project 2022/Estimation/Data/publication_calendar_final_final.xlsx','Sheet','Variables');
macro_groups = T2.Group(1:size(output.incl,2));
macro_groups(17) = 6; % Changed the VISA grouping manually for now


%% Inclusion Probability Graphs
startd = find(d_q == '31-Dec-2019');
pincl_almon = almon.incl;
pincl_hst_sw = hs_sw_svt.incl;
pincl_hs_sw = hs_sw_sv.incl;
% Posterior Inclusion Probability Stacked Bar Chart (Average over all time periods)
dat_almon = squeeze(mean(pincl_almon(:,:,1:startd),3));
dat_hst_sw = squeeze(mean(pincl_hst_sw(:,:,1:startd-1),3));
dat_hs_sw = squeeze(mean(pincl_hs_sw(:,:,1:startd-1),3));
%dat_hs(20:end,9) = dat_hs(20:end,9)-0.15; 


%% Create Combined Heatmap: Pre Pandemic Average

fig1 = figure;

YLabels = 1:20;
% Convert each number in the array into a string
CustomYLabels = string(YLabels);
% Replace all but the fifth elements by spaces
CustomYLabels(mod(YLabels,3) ~= 0) = " ";
% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels



subplot(3,1,1)
hm = heatmap(dat_almon(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none')
colormap(flipud(hot))
hm.GridVisible = 'off';
hm.XDisplayLabels = {varnames_qm{1:17}};
hm.YDisplayLabels = CustomYLabels;
ylabel('Nowcast Periods')
title('KP (2023)')
hm.FontSize = 10;


subplot(3,1,2)
hm = heatmap(dat_hst_sw(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none')
colormap(flipud(hot))
hm.GridVisible = 'off';
hm.XDisplayLabels = {varnames_qm{1:17}};
hm.YDisplayLabels = CustomYLabels;
ylabel('Nowcast Periods')
title('Trend(HS,no g)-SV(HS)-t-GIGG')
hm.FontSize = 10;

subplot(3,1,3)
hm = heatmap(dat_hs_sw(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none')
colormap(flipud(hot))
hm.GridVisible = 'off';
hm.XDisplayLabels = {varnames_qm{1:17}};
hm.YDisplayLabels = CustomYLabels;
ylabel('Nowcast Periods')
title('KP (2023) - hierarchical')
hm.FontSize = 10;

set(fig1,'Position',get(0,'Screensize'))
h=fig1;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

sgtitle('Average Inclusion Probability Pre-Pandemic',FontSize = 16);
extraname = 'pre_pandemic_newmods_';
saveas(h,strcat(savedir,extraname,modelname1))




%% Combined Heatmap: Post Pandemic Average
startd = find(d_q == '31-Dec-2019');
pincl_almon = almon.incl;
pincl_hst_sw = hs_sw_svt.incl;
% Posterior Inclusion Probability Stacked Bar Chart (Average over all time periods)
dat_almon = squeeze(mean(pincl_almon(:,:,startd+1:end),3));
dat_hst_sw = squeeze(mean(pincl_hst_sw(:,:,startd:end),3));
dat_hs_sw = squeeze(mean(pincl_hs_sw(:,:,startd:end),3));
%dat_hs(3:end,9) = dat_hs(3:end,9)-0.15; 


fig2 = figure;

YLabels = 1:20;
% Convert each number in the array into a string
CustomYLabels = string(YLabels);
% Replace all but the fifth elements by spaces
CustomYLabels(mod(YLabels,3) ~= 0) = " ";
% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels



subplot(3,1,1)
hm = heatmap(dat_almon(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none')
colormap(flipud(hot))
hm.GridVisible = 'off';
hm.XDisplayLabels = {varnames_qm{1:17}};
hm.YDisplayLabels = CustomYLabels;
ylabel('Nowcast Periods')
title('KP (2023)')
hm.FontSize = 10;


subplot(3,1,2)
hm = heatmap(dat_hst_sw(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none')
colormap(flipud(hot))
hm.GridVisible = 'off';
hm.XDisplayLabels = {varnames_qm{1:17}};
hm.YDisplayLabels = CustomYLabels;
ylabel('Nowcast Periods')
title('Trend(HS)-SV(HS)-t-GIGG')
hm.FontSize = 10;

subplot(3,1,3)
hm = heatmap(dat_hs_sw(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none')
colormap(flipud(hot))
hm.GridVisible = 'off';
hm.XDisplayLabels = {varnames_qm{1:17}};
hm.YDisplayLabels = CustomYLabels;
ylabel('Nowcast Periods')
title('Trend(HS)-SV(HS)-GIGG')
hm.FontSize = 10;

set(fig2,'Position',get(0,'Screensize'))
h=fig2;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

sgtitle('Average Inclusion Probability Including the Pandemic',FontSize = 16);
extraname = 'post_pandemic_newmods_';
saveas(h,strcat(savedir,extraname,modelname1))







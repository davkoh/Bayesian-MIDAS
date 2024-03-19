%% In-Sample Result Figures: V2 with larger confidence sets %%
    % Trend
    % SV-trend
    % SV observation equation
    % Cyclical component
clear;
clear all

%% Upload the Data
cd '/Users/dk/Library/CloudStorage/OneDrive-Heriot-WattUniversity/Education/BoE Nowcasting Project 2022/Scripts from PC/Downloads/GIGG-Nowcast-Output' 
addpath('/Users/dk/Library/CloudStorage/OneDrive-Heriot-WattUniversity/Education/BoE Nowcasting Project 2022/Functions')
addpath('/Users/dk/Downloads')

% SV's
load('trendsv2.mat')
SV_trend = exp(1/2*out.g);
SV_obs = exp(1/2*out.h);

% Trend & Cyclical Component
trend = out.tau;
cyc = out.cyc;

% GDP series
load('matlab.mat')
y = yin(1:end-1);

%% V1: Pre & Post Pandemic trends in the upper panels, SV & Trends in lower panels

% Upload a results matrix for the dates needed
load(['/Users/dk/Library/CloudStorage/OneDrive-Heriot-WattUniversity/Education' ...
    '/BoE Nowcasting Project 2022/Estimation/Data' ...
    '/UKdata_11_11_1980_2021_GDPCONINVHOUR_sur_act_lab_mort_ie_vis.mat']);
nfor = 42;
dq = d_q(end-size(y)+1:end); 

% Adjust the series
SVobs_mean = mean(SV_obs,2);
SVobs_lower = quantile(SV_obs',0.05);
SVobs_upper = quantile(SV_obs',0.95);

SVtrend_mean = mean(SV_trend,2);
SVtrend_lower = quantile(SV_trend',0.05);
SVtrend_upper = quantile(SV_trend',0.95);

trend_mean = mean(trend,2); %mean(trend,2)-mean(mean(trend,2))+ mean(y);
trend_lower = quantile(trend',0.12);
trend_upper = quantile(trend',0.88);

cyc_mean = mean(cyc,2);%-mean(mean(cyc,2)) + mean(y); %mean(trend,2)-mean(mean(trend,2))+ mean(y);
cyc_lower = quantile(cyc',0.05);% - mean(quantile(cyc',0.2)) + mean(y);
cyc_upper = quantile(cyc',0.95);% - mean(mean(cyc,2)) + mean(y);


fig1 = figure;

startdate = datenum(dq(1));
enddate_pre = datenum(dq(end-5)); 
startdate_post = datenum(dq(end-5));
enddate_post = datenum(dq(end));
xData_pre = datenum(dq(1:end-5))';
xData_post = datenum(dq(end-5:end))';

hAx(1) = subplot(2,2,1);
plot(xData_pre,cyc_mean(1:end-5),'LineWidth',1.5)
hold on
[ph,msg]=jbfill(xData_pre,cyc_upper((1:end-5)),cyc_lower((1:end-5)),[0 0.4470 0.7410],'none',1,0.3);
hold on
plot(xData_pre,trend_mean(1:end-5),'LineWidth',1.5)
[ph,msg]=jbfill(xData_pre,trend_upper((1:end-5)),trend_lower((1:end-5)),[0.8500 0.3250 0.0980],'none',1,0.3);
hold on
plot(xData_pre,y(1:end-5,:),'LineWidth',1.5,'Color','black','LineStyle','-.')
datetick('x','QQ-YY','keepticks')
xlim([startdate enddate_pre]);
xticks(xData_pre(3:6:end))
xticklabels(datestr(xData_pre(3:6:end),'QQ-YY'))
title('Trend & Cyclical: Pre-Pandemic',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(1),'XTickLabelRotation',45)

hAx(2) = subplot(2,2,2);
plot(xData_post,cyc_mean(end-5:end),'LineWidth',1.5)
hold on
[ph,msg]=jbfill(xData_post,cyc_upper(end-5:end),cyc_lower((end-5:end)),[0 0.4470 0.7410],'none',1,0.3);
hold on
plot(xData_post,trend_mean(end-5:end),'LineWidth',1.5)
[ph,msg]=jbfill(xData_post,trend_upper((end-5:end)),trend_lower((end-5:end)),[0.8500 0.3250 0.0980],'none',1,0.3);
hold on
plot(xData_post,y(end-5:end,:),'LineWidth',1.5,'Color','black','LineStyle','-.')
datetick('x','QQ-YY','keepticks')
xlim([startdate_post enddate_post]);
xticks(xData_post)
xticklabels(datestr(xData_post,'QQ-YY'))
title('Trend & Cyclical: Pandemic',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(2),'XTickLabelRotation',45)

hAx(3) = subplot(2,2,3);
plot(datenum(dq)',SVobs_mean,'LineWidth',1.5)
hold on
[ph,msg]=jbfill(datenum(dq)',SVobs_upper,SVobs_lower,[0 0.4470 0.7410],'none',1,0.3);
hold on
datetick('x','QQ-YY','keepticks')
xlim([startdate enddate_post]);
xticks(datenum(dq(3:6:end)))
xticklabels(datestr(datenum(dq(3:6:end)),'QQ-YY'))
title('SV: Observation',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(3),'XTickLabelRotation',45)

hAx(4) = subplot(2,2,4);
plot(datenum(dq)',SVtrend_mean,'LineWidth',1.5)
hold on
[ph,msg]=jbfill(datenum(dq)',SVtrend_upper,SVtrend_lower,[0 0.4470 0.7410],'none',1,0.3);
hold on
datetick('x','QQ-YY','keepticks')
xlim([startdate enddate_post]);
ylim([0 0.2])
xticks(datenum(dq(3:6:end)))
xticklabels(datestr(datenum(dq(3:6:end)),'QQ-YY'))
title('SV: Trend',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(4),'XTickLabelRotation',45)

%hL =  legend('','U-BMIDAS Combination','Trend-BMIDAS-HS','M&S','Trend-BMIDAS-GIGG','Realised');
% Deal with Legend & Readability
%newPosition = [0.05,0.8,0.03,0.098406745154051];
% Because I would need two different legends, I have not inlcuded one for
% now
%set(hL,'Position', newPosition,'FontSize',8);
h=fig1;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

figname = ['/Users/dk/Documents/GitHub/Bayesian-MIDAS/Output/in_sample_trends.pdf'];
saveas(fig1,figname)


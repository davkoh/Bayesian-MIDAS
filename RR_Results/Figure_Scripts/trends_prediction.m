%% In-Sample Result Figures: V2 with larger confidence sets %%
    % Trend
    % SV-trend
    % SV observation equation
    % Cyclical component
clear;
clear all



load("../Data/UK_dat_2024.mat");

mod1 = load("Output_iterated/results_iteratednowcasts_trendssave_newgigg_oldsv_bg_0.5_trend_sv_almon_ortho_groupsparse.mat");
mod1 = mod1.output;

periods = [6 18 26];

%% Get the trends and their respective quantiles
% Trend
clear trend
for i = 1:length(periods)
    trend(i).trend = squeeze(mod1.tau_all(:,periods(i),:)) - mean(squeeze(median(mod1.tau_all(:,periods(i),:),1))) + mean(y);
    trend(i).loc = median(trend(i).trend,1); 
    trend(i).trend_lower= quantile(trend(i).trend,0.05); 
    trend(i).trend_upper= quantile(trend(i).trend,0.95); 
end

% SV
clear SV
for i = 1:length(periods)
    SV(i).SV = exp(1/2*squeeze(mod1.sv_all(:,periods(i),:)));
    SV(i).loc = median(SV(i).SV,1); 
    SV(i).SV_lower= quantile(SV(i).SV,0.05); 
    SV(i).SV_upper= quantile(SV(i).SV,0.95); 
end


% SV of Trend
clear SV_Trend
for i = 1:length(periods)
    SV_Trend(i).SV_Trend = exp(1/2*squeeze(mod1.sv_trend_all(:,periods(i),:)));
    SV_Trend(i).loc = median(SV_Trend(i).SV_Trend,1); 
    SV_Trend(i).trend_lower= quantile(SV_Trend(i).SV_Trend,0.05); 
    SV_Trend(i).trend_upper= quantile(SV_Trend(i).SV_Trend,0.95); 
end

% Cyclical
clear Cyclical
for i = 1:length(periods)
    Cyclical(i).Cyclical = squeeze(mod1.cyc_pred_all(:,periods(i),:)) - mean(squeeze(median(mod1.cyc_pred_all(:,periods(i),:),1))) + mean(y);
    Cyclical(i).loc = median(Cyclical(i).Cyclical,1); 
    Cyclical(i).Cyclical_lower= quantile(Cyclical(i).Cyclical,0.05); 
    Cyclical(i).Cyclical_upper= quantile(Cyclical(i).Cyclical,0.95); 
end

%% Plot the Figure

dq = d_q;
dq = dq(end-nfor+1:end);
y = y(end-nfor+1:end);

fig1 = figure;

startdate = datenum(dq(1));
enddate_pre = datenum(dq(end-13)); 
startdate_post = datenum(dq(end-12));
enddate_post = datenum(dq(end));
xData_pre = datenum(dq(1:end-13))';
xData_post = datenum(dq(end-12:end))';



hAx(2) = subplot(3,2,2);
plot(xData_post,Cyclical(1).loc(end-12:end),'LineWidth',1.5,Color="#8F968B",LineStyle="-")
hold on
plot(xData_post,Cyclical(1).Cyclical_lower(end-12:end),'LineWidth',1,Color="#8F968B",LineStyle=":")
plot(xData_post,Cyclical(1).Cyclical_upper(end-12:end),'LineWidth',1,Color="#8F968B",LineStyle=":")

plot(xData_post,Cyclical(2).loc(end-12:end),'LineWidth',1.5,Color="#3EA772",LineStyle="-")
plot(xData_post,Cyclical(2).Cyclical_lower(end-12:end),'LineWidth',1,Color="#3EA772",LineStyle=":")
plot(xData_post,Cyclical(2).Cyclical_upper(end-12:end),'LineWidth',1,Color="#3EA772",LineStyle=":")


plot(xData_post,y(end-12:end,:),'LineWidth',1,'Color','blue','LineStyle','-')

plot(xData_post,Cyclical(3).loc(end-12:end),'LineWidth',1.5,Color="#000000",LineStyle="-")
plot(xData_post,Cyclical(3).Cyclical_lower(end-12:end),'LineWidth',1,Color="#000000",LineStyle=":")
plot(xData_post,Cyclical(3).Cyclical_upper(end-12:end),'LineWidth',1,Color="#000000",LineStyle=":")

datetick('x','QQ-YY','keepticks')
xlim([startdate_post enddate_post]);
xticks(xData_post)
xticklabels(datestr(xData_post,'QQ-YY'))
title('Cycle: Pandemic',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(2),'XTickLabelRotation',45)

hAx(3) = subplot(3,2,3);

plot(xData_pre,trend(1).loc(1:end-13),'LineWidth',1.5,Color="#8F968B",LineStyle="-")
hold on
plot(xData_pre,trend(1).trend_lower(1:end-13),'LineWidth',1,Color="#8F968B",LineStyle=":")
plot(xData_pre,trend(1).trend_upper(1:end-13),'LineWidth',1,Color="#8F968B",LineStyle=":")

plot(xData_pre,trend(2).loc(1:end-13),'LineWidth',1.5,Color="#3EA772",LineStyle="-")
plot(xData_pre,trend(2).trend_lower(1:end-13),'LineWidth',1,Color="#3EA772",LineStyle=":")
plot(xData_pre,trend(2).trend_upper(1:end-13),'LineWidth',1,Color="#3EA772",LineStyle=":")


plot(xData_pre,trend(3).loc(1:end-13),'LineWidth',1.5,Color="#000000",LineStyle="-")
plot(xData_pre,trend(3).trend_lower(1:end-13),'LineWidth',1,Color="#000000",LineStyle=":")
plot(xData_pre,trend(3).trend_upper(1:end-13),'LineWidth',1,Color="#000000",LineStyle=":")

datetick('x','QQ-YY','keepticks')
%ylim([-1 2])
xlim([startdate enddate_pre]);
xticks(xData_pre(3:6:end))
xticklabels(datestr(xData_pre(3:6:end),'QQ-YY'))
title('Trend: Pre-Pandemic',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(3),'XTickLabelRotation',45)

hAx(4) = subplot(3,2,4);


plot(xData_post,trend(1).loc(end-12:end),'LineWidth',1.5,Color="#8F968B",LineStyle="-")
hold on
plot(xData_post,trend(1).trend_lower(end-12:end),'LineWidth',1,Color="#8F968B",LineStyle=":")
plot(xData_post,trend(1).trend_upper(end-12:end),'LineWidth',1,Color="#8F968B",LineStyle=":")

plot(xData_post,trend(2).loc(end-12:end),'LineWidth',1.5,Color="#3EA772",LineStyle="-")
plot(xData_post,trend(2).trend_lower(end-12:end),'LineWidth',1,Color="#3EA772",LineStyle=":")
plot(xData_post,trend(2).trend_upper(end-12:end),'LineWidth',1,Color="#3EA772",LineStyle=":")


plot(xData_post,trend(3).loc(end-12:end),'LineWidth',1.5,Color="#000000",LineStyle="-")
plot(xData_post,trend(3).trend_lower(end-12:end),'LineWidth',1,Color="#000000",LineStyle=":")
plot(xData_post,trend(3).trend_upper(end-12:end),'LineWidth',1,Color="#000000",LineStyle=":")


datetick('x','QQ-YY','keepticks')
xlim([startdate_post enddate_post]);
xticks(xData_post)
ylim([-2 2])
xticklabels(datestr(xData_post,'QQ-YY'))
title('Trend: Pandemic',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(4),'XTickLabelRotation',45)

hAx(5) = subplot(3,2,5);

plot(datenum(dq)',SV(1).loc,'LineWidth',1.5,Color="#8F968B",LineStyle="-")
hold on
plot(datenum(dq)',SV(1).SV_lower,'LineWidth',1,Color="#8F968B",LineStyle=":")
plot(datenum(dq)',SV(1).SV_upper,'LineWidth',1,Color="#8F968B",LineStyle=":")

plot(datenum(dq)',SV(2).loc,'LineWidth',1.5,Color="#3EA772",LineStyle="-")
plot(datenum(dq)',SV(2).SV_lower,'LineWidth',1,Color="#3EA772",LineStyle=":")
plot(datenum(dq)',SV(2).SV_upper,'LineWidth',1,Color="#3EA772",LineStyle=":")


plot(datenum(dq)',SV(3).loc,'LineWidth',1.5,Color="#000000",LineStyle="-")
plot(datenum(dq)',SV(3).SV_lower,'LineWidth',1,Color="#000000",LineStyle=":")
plot(datenum(dq)',SV(3).SV_upper,'LineWidth',1,Color="#000000",LineStyle=":")

datetick('x','QQ-YY','keepticks')
xlim([startdate enddate_post]);
ylim([0 15])
xticks(datenum(dq(3:6:end)))
xticklabels(datestr(datenum(dq(3:6:end)),'QQ-YY'))
title('SV: Observation',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(5),'XTickLabelRotation',45)


hAx(6) = subplot(3,2,6);
plot(datenum(dq)',SV_Trend(1).loc,'LineWidth',1.5,Color="#8F968B",LineStyle="-")
hold on
plot(datenum(dq)',SV_Trend(1).trend_lower,'LineWidth',1,Color="#8F968B",LineStyle=":")
plot(datenum(dq)',SV_Trend(1).trend_upper,'LineWidth',1,Color="#8F968B",LineStyle=":")

plot(datenum(dq)',SV_Trend(2).loc,'LineWidth',1.5,Color="#3EA772",LineStyle="-")
plot(datenum(dq)',SV_Trend(2).trend_lower,'LineWidth',1,Color="#3EA772",LineStyle=":")
plot(datenum(dq)',SV_Trend(2).trend_upper,'LineWidth',1,Color="#3EA772",LineStyle=":")


plot(datenum(dq)',SV_Trend(3).loc,'LineWidth',1.5,Color="#000000",LineStyle="-")
plot(datenum(dq)',SV_Trend(3).trend_lower,'LineWidth',1,Color="#000000",LineStyle=":")
plot(datenum(dq)',SV_Trend(3).trend_upper,'LineWidth',1,Color="#000000",LineStyle=":")

datetick('x','QQ-YY','keepticks')
xlim([startdate enddate_post]);
%ylim([0 1])
xticks(datenum(dq(3:6:end)))
xticklabels(datestr(datenum(dq(3:6:end)),'QQ-YY'))
title('SV: Trend',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(6),'XTickLabelRotation',45)

hAx(1) = subplot(3,2,1);

plot(xData_pre,Cyclical(1).loc(1:end-13),'LineWidth',1.5,Color="#8F968B",LineStyle="-")
hold on
plot(xData_pre,Cyclical(1).Cyclical_lower(1:end-13),'LineWidth',1,Color="#8F968B",LineStyle=":")
plot(xData_pre,Cyclical(1).Cyclical_upper(1:end-13),'LineWidth',1,Color="#8F968B",LineStyle=":")

plot(xData_pre,Cyclical(2).loc(1:end-13),'LineWidth',1.5,Color="#3EA772",LineStyle="-")
plot(xData_pre,Cyclical(2).Cyclical_lower(1:end-13),'LineWidth',1,Color="#3EA772",LineStyle=":")
plot(xData_pre,Cyclical(2).Cyclical_upper(1:end-13),'LineWidth',1,Color="#3EA772",LineStyle=":")

plot(xData_pre,Cyclical(3).loc(1:end-13),'LineWidth',1.5,Color="#000000",LineStyle="-")
plot(xData_pre,Cyclical(3).Cyclical_lower(1:end-13),'LineWidth',1,Color="#000000",LineStyle=":")
plot(xData_pre,Cyclical(3).Cyclical_upper(1:end-13),'LineWidth',1,Color="#000000",LineStyle=":")


plot(xData_pre,y(1:end-13,:),'LineWidth',1,'Color','blue','LineStyle','-')

datetick('x','QQ-YY','keepticks')
xlim([startdate enddate_pre]);
xticks(xData_pre(3:6:end))
xticklabels(datestr(xData_pre(3:6:end),'QQ-YY'))
title('Cycle: Pre-Pandemic',FontSize=16)
h=get(fig1,'CurrentAxes')
set(h,'FontSize',15)
set(hAx(1),'XTickLabelRotation',45)

hL =  legend('Period:1','','','Period:12','','','Period:19','','','GDP');

newPosition = [0.05,0.8,0.03,0.098406745154051];
% Because I would need two different legends, I have not inlcuded one for
% now
set(hL,'Position', newPosition,'FontSize',10);
h=fig1;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);




modname = strcat('Figures/insampletrends','_newgigg_trend_sv_predictions_','all',".pdf");

saveas(fig1,modname)


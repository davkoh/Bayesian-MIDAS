%% Inclusion probability Map

% Define the model directory
load('Output/TREND_PC_OBS_PC_PRIOR_gigg_GIGGTYPE_fixed_TAIL_normal_TRANSFORM_almon/results.mat');

% Automatic from here
mod = output;
var_names = output.var_names;
dq_nfor = output.d_q;
dqend = output.d_q(end);
pincl = output.incl;

figname = mod.mod_name;

%% Inclusion Probability Graph
% Start Dates of sub-samples
startDate = datetime('01-Jan-2007');
endDate = datetime('31-Dec-2009');
dq_dates = dateshift(dq_nfor, 'start', 'day');
indices_gfc = find(dq_dates >= startDate & dq_dates <= endDate)';

startDate = datetime('01-Jan-2010');
endDate = datetime('31-Dec-2019');
indices_tranq = find(dq_dates >= startDate & dq_dates <= endDate)';

startDate = datetime('01-Jan-2020');
endDate = dqend;
indices_covid = find(dq_dates >= startDate & dq_dates <= endDate)';

% Retrieve inclusion probabilities
dat_fullsamp = squeeze(mean(pincl(:,:,:),3));
dat_GFC= squeeze(mean(pincl(:,:,indices_gfc),3));
dat_Tranq= squeeze(mean(pincl(:,:,indices_tranq),3));
dat_Covid= squeeze(mean(pincl(:,:,indices_covid),3));

%% Heatmap: Periods


fig1 = figure;

YLabels = 1:19;
% Convert each number in the array into a string
CustomYLabels = string(YLabels);
CustomYLabels([2 3 5 6 8 9 11 12 14 15 17 18]) = "";

subplot(3,1,1)
hm = heatmap(dat_fullsamp(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none')
colormap(flipud(hot))
hm.GridVisible = 'off';
hm.XDisplayLabels = var_names;
hm.YDisplayLabels = CustomYLabels;
ylabel('Nowcast Periods')
title('Full sample (2007q3-2023q3)')
hm.FontSize = 15;

subplot(3,1,2)
hm = heatmap(dat_GFC(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none')
colormap(flipud(hot))
hm.GridVisible = 'off';
hm.XDisplayLabels = var_names;
hm.YDisplayLabels = CustomYLabels;
ylabel('Nowcast Periods')
title('Pre-pandemic (2007q3-2019q4)')
hm.FontSize = 15;


subplot(3,1,3)
hm = heatmap(dat_Covid(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none')
colormap(flipud(hot))
hm.GridVisible = 'off';
hm.XDisplayLabels = var_names;
hm.YDisplayLabels = CustomYLabels;
ylabel('Nowcast Periods')
title('Pandemic (2020q1-2023q3)')
hm.FontSize = 15;

set(fig1,'Position',get(0,'Screensize'))
h=fig1;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [-0.05 0 1.1 1]);

saveas(h, fullfile("Figures", figname + ".fig"));
saveas(h, fullfile("Figures", figname + ".pdf"));






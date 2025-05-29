function plot_heatmap_function(resultsFile)
% plotInclusionProbabilityMap - Plots inclusion probability heatmaps for different periods.
%
% Usage:
%   plotInclusionProbabilityMap('Output/TREND_PC_OBS_PC_PRIOR_gigg_GIGGTYPE_fixed_TAIL_normal_TRANSFORM_almon/results.mat')
%
% Inputs:
%   resultsFile - Path to the .mat file containing the results (string)

figuresDir = 'Figures';

    % Load results
    load(resultsFile, 'output');

    mod = output;
    var_names = output.var_names;
    dq_nfor = output.d_q;
    dqend = output.d_q(end);
    pincl = output.incl;

    figname = mod.mod_name;

    % Dates for sub-samples
    dq_dates = dateshift(dq_nfor, 'start', 'day');
    
    % Pre covid period
    indices_pre = find(dq_dates >= datetime('01-Jan-2007') & dq_dates <= datetime('31-Dec-2019'))';
    % GFC period
    indices_gfc = find(dq_dates >= datetime('01-Jan-2007') & dq_dates <= datetime('31-Dec-2009'))';
    % Tranquil period
    indices_tranq = find(dq_dates >= datetime('01-Jan-2010') & dq_dates <= datetime('31-Dec-2019'))';
    % Covid period
    indices_covid = find(dq_dates >= datetime('01-Jan-2020') & dq_dates <= dqend)';

    % Retrieve inclusion probabilities
    dat_fullsamp = squeeze(mean(pincl(:,:,:),3));
    dat_pre = squeeze(mean(pincl(:,:,indices_pre),3));
    dat_GFC = squeeze(mean(pincl(:,:,indices_gfc),3));
    dat_Tranq = squeeze(mean(pincl(:,:,indices_tranq),3));
    dat_Covid = squeeze(mean(pincl(:,:,indices_covid),3));

    % Heatmap: Periods
    fig1 = figure;

    YLabels = 1:19;
    CustomYLabels = string(YLabels);
    CustomYLabels([2 3 5 6 8 9 11 12 14 15 17 18]) = "";

    subplot(2,1,1)
    hm = heatmap(dat_fullsamp(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none');
    colormap(flipud(hot))
    hm.GridVisible = 'off';
    hm.XDisplayLabels = var_names;
    hm.YDisplayLabels = CustomYLabels;
    ylabel('Nowcast Periods')
    title('Pre-pandemic (2007q3-2019q4)')
    hm.FontSize = 15;

    %{
    subplot(3,1,2)
    hm = heatmap(dat_GFC(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none');
    colormap(flipud(hot))
    hm.GridVisible = 'off';
    hm.XDisplayLabels = var_names;
    hm.YDisplayLabels = CustomYLabels;
    ylabel('Nowcast Periods')
    title('Pre-pandemic (2007q3-2019q4)')
    hm.FontSize = 15;
    %}

    subplot(2,1,2)
    hm = heatmap(dat_Covid(1:end,:),'ColorLimits',[0 1],'CellLabelColor','none');
    colormap(flipud(hot))
    hm.GridVisible = 'off';
    hm.XDisplayLabels = var_names;
    hm.YDisplayLabels = CustomYLabels;
    ylabel('Nowcast Periods')
    title('Pandemic (2020q1-2023q3)')
    hm.FontSize = 15;

    set(fig1,'Position',get(0,'Screensize'))
    set(fig1,'PaperOrientation','landscape');
    set(fig1,'PaperUnits','normalized');
    set(fig1,'PaperPosition', [-0.05 0 1.1 1]);

    if ~exist(figuresDir, 'dir')
        mkdir(figuresDir);
    end

    saveas(fig1, fullfile(figuresDir, figname + ".fig"));
    %saveas(fig1, fullfile(figuresDir, figname + ".pdf"));
end

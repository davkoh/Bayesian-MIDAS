function plot_trend_decomposition(output_file, periods,foldername)

% inputs: ouput including estimates for trend and SVs in trend and
% observation equation
% nowcast periods (between 1 and vint) to plot posterior estimates for
% foldername: outputfolder for figure

% Load model data
out=output_file;

% Extract necessary data
y = out.yf;
dq = out.d_q;

% Define pandemic cutoff date
pandemic_cutoff_date = datenum('01-Jan-2020');

% Determine pre- and post-pandemic indices
pre_pandemic_indices = find(datenum(dq) < pandemic_cutoff_date);
post_pandemic_indices = find(datenum(dq) >= pandemic_cutoff_date);

% Date ranges for plotting
xData_pre = datenum(dq(pre_pandemic_indices))';
xData_post = datenum(dq(post_pandemic_indices))';
startdate = min(xData_pre);
enddate_pre = max(xData_pre);
enddate_post = max(datenum(dq));

% Process trends and their respective quantiles
trends = get_trends(out, periods, y);
SV = get_SV(out, periods);
SV_Trend = get_SV_Trend(out, periods);
cyclical = get_cyclical(out, periods, y);

% Create Figure
fig = figure;


create_subplot(fig, xData_post, cyclical, y, 'Cycle: Pandemic', [3, 2, 2], post_pandemic_indices, [min(xData_post), max(xData_post)], xData_post, true, [-30, 22], 'yes');
create_subplot(fig, xData_pre, trends, y, 'Trend: Pre-Pandemic', [3, 2, 3], pre_pandemic_indices, [min(xData_pre), max(xData_pre)], xData_pre(3:6:end), true, [-2, 2], 'no');
create_subplot(fig, xData_post, trends, y, 'Trend: Pandemic', [3, 2, 4], post_pandemic_indices, [min(xData_post), max(xData_post)], xData_post, false, [], 'no');
create_subplot(fig, datenum(dq)', SV, y, 'SV: Observation', [3, 2, 5], 1:length(dq), [startdate, enddate_post], datenum(dq(3:6:end)), true, [0, 15], 'no');
create_subplot(fig, datenum(dq)', SV_Trend, y, 'SV: Trend', [3, 2, 6], 1:length(dq), [startdate, enddate_post], datenum(dq(3:6:end)), false, [], 'no');
create_subplot(fig, xData_pre, cyclical, y, 'Cycle: Pre-Pandemic', [3, 2, 1], pre_pandemic_indices, [min(xData_pre), max(xData_pre)], xData_pre(3:6:end), true, [-2.5, 2], 'yes');

% Add a legend to the first plot (Cycle: Pre-Pandemic)
legend_labels = {'Period: 1', '', '', 'Period: 12', '', '', 'Period: 19', '','', 'GDP'};
lgd = legend(legend_labels, 'Orientation', 'vertical', 'NumColumns', 1,FontSize=8);
set(lgd, 'Position', [0.001, 0.75, 0.1, 0.1]); % Adjust the position as needed to ensure no touching

% Set figure properties for saving
set(fig, 'PaperOrientation', 'landscape');
set(fig, 'PaperUnits', 'normalized');
set(fig, 'PaperPosition', [0 0 1 1]);


% Save figure
modname = fullfile(foldername, strcat('Figure_Trend_decomp.pdf'));
exportgraphics(fig, modname, 'ContentType', 'vector', 'BackgroundColor', 'white', 'Resolution', 300);
saveas(fig, strrep(modname, '.pdf', '.fig'));

end

function trends = get_trends(mod, periods, y)
    % Extract trend data from the model
    for i = 1:length(periods)
        t = squeeze(mod.tau_all(:, periods(i), :));
        trends(i).trend = t - mean(squeeze(median(t, 1))) + mean(y);
        trends(i).loc = median(trends(i).trend, 1); 
        trends(i).lower = quantile(trends(i).trend, 0.05); 
        trends(i).upper = quantile(trends(i).trend, 0.95);
    end
end

function SV = get_SV(mod, periods)
    % Extract SV data from the model
    for i = 1:length(periods)
        s = exp(1/2 * squeeze(mod.sv_all(:, periods(i), :)));
        SV(i).SV = s;
        SV(i).loc = median(SV(i).SV, 1); 
        SV(i).lower = quantile(SV(i).SV, 0.05); 
        SV(i).upper = quantile(SV(i).SV, 0.95);
    end
end

function SV_Trend = get_SV_Trend(mod, periods)
    % Extract SV Trend data from the model
    for i = 1:length(periods)
        st = exp(1/2 * squeeze(mod.sv_trend_all(:, periods(i), :)));
        SV_Trend(i).SV_Trend = st;
        SV_Trend(i).loc = median(SV_Trend(i).SV_Trend, 1); 
        SV_Trend(i).lower = quantile(SV_Trend(i).SV_Trend, 0.05); 
        SV_Trend(i).upper = quantile(SV_Trend(i).SV_Trend, 0.95); 
    end
end

function cyclical = get_cyclical(mod, periods, y)
    % Extract cyclical data from the model
    for i = 1:length(periods)
        c = squeeze(mod.cyc_pred_all(:, periods(i), :));
        cyclical(i).Cyclical = c - mean(squeeze(median(c, 1))) + mean(y);
        cyclical(i).loc = median(cyclical(i).Cyclical, 1); 
        cyclical(i).lower = quantile(cyclical(i).Cyclical, 0.05); 
        cyclical(i).upper = quantile(cyclical(i).Cyclical, 0.95);
    end
end

function create_subplot(fig, xData, structure, y, titleText, subplot_pos, indices, xlim_range, xticks_range, use_ylim, ylim_range, include_y)
    % Create subplot with confidence intervals and median predictions
    hAx = subplot(subplot_pos(1), subplot_pos(2), subplot_pos(3), 'Parent', fig);
    hold(hAx, 'on')
    periodColors = ["#8F968B", "#3EA772", "#000000"];
    
    for i = 1:length(structure)
        plot(hAx, xData, structure(i).loc(indices), 'LineWidth', 1.5, 'Color', periodColors(i));
        plot(hAx, xData, structure(i).lower(indices), 'LineWidth', 1, 'Color', periodColors(i), 'LineStyle', ':');
        plot(hAx, xData, structure(i).upper(indices), 'LineWidth', 1, 'Color', periodColors(i), 'LineStyle', ':');
    end
    
    if strcmp(include_y, 'yes')
        plot(hAx, xData, y(indices,:), 'LineWidth', 1, 'Color', 'blue', 'LineStyle', '-');
    end
    
    datetick(hAx, 'x', 'QQ-YY', 'keepticks');
    xlim(hAx, xlim_range);
    xticks(hAx, xticks_range);
    xticklabels(hAx, datestr(xticks_range, 'QQ-YY'))
    
    if use_ylim && ~isempty(ylim_range)
        ylim(hAx, ylim_range);
    end
    
    title(hAx, titleText, 'FontSize', 12);
    set(hAx, 'FontSize', 10);
    set(hAx, 'XTickLabelRotation', 45);
    hold(hAx, 'off');
end

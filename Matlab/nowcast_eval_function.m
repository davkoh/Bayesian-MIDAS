function nowcast_eval_function(model_files)
%NOWCAST_EVAL Plots RMSFE and CRPS prediction graphs for multiple models.
%   nowcast_eval(model_files)
%   model_files: cell array of strings, each is a path to a model .mat file

% Automatically extract model names from each model's output
model_names = cell(1, numel(model_files));
for m = 1:numel(model_files)
    tmp = load(model_files{m});
    if isfield(tmp.output, 'mod_name')
        model_names{m} = tmp.output.mod_name;
    else
        model_names{m} = ['Model ' num2str(m)];
    end
end
nModels = numel(model_files);

% Preallocate cell arrays for results
mod_output = cell(1, nModels);
rmsfe_gfc = cell(1, nModels);
rmsfe_tranq = cell(1, nModels);
rmsfe_covid = cell(1, nModels);
rtcrps_gfc = cell(1, nModels);
rtcrps_tranq = cell(1, nModels);
rtcrps_covid = cell(1, nModels);

% Load models and compute metrics
for m = 1:nModels
    tmp = load(model_files{m});
    mod_output{m} = tmp.output;
    res_mod = mod_output{m}.resid_all;
    crps_mod = mod_output{m}.crps_all;
    % Use the same dq_nfor and dqend for all models
    if m == 1
        dq_nfor = mod_output{m}.d_q;
        dqend = mod_output{m}.d_q(end);
    end

    % GFC
    startDate = datetime('01-Jan-2007');
    endDate = datetime('31-Dec-2009');
    dq_dates = dateshift(dq_nfor, 'start', 'day');
    indices_gfc = find(dq_dates >= startDate & dq_dates <= endDate)';
    rmsfe_gfc{m} = std(res_mod(:,indices_gfc)')';
    rtcrps_gfc{m} = mean(crps_mod(:,indices_gfc),2);

    % Tranquil
    startDate = datetime('01-Jan-2010');
    endDate = datetime('31-Dec-2019');
    indices_tranq = find(dq_dates >= startDate & dq_dates <= endDate)';
    rmsfe_tranq{m} = std(res_mod(:,indices_tranq)')';
    rtcrps_tranq{m} = mean(crps_mod(:,indices_tranq),2);

    % Covid
    startDate = datetime('01-Jan-2020');
    endDate = dqend;
    indices_covid = find(dq_dates >= startDate & dq_dates <= endDate)';
    rmsfe_covid{m} = std(res_mod(:,indices_covid)')';
    rtcrps_covid{m} = mean(crps_mod(:,indices_covid),2);
end

%% Figure: Eval Graph for all models

fig1 = figure;
t = tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact'); % Use tiledlayout

colors = lines(nModels); % Distinct colors for each model

% First row: RMSFE
for s = 1:3
    ax = nexttile(s);
    hold(ax, 'on')
    for m = 1:nModels
        switch s
            case 1
                plot(ax, rmsfe_gfc{m}, 'LineWidth',2, 'Color',colors(m,:), 'LineStyle','-');
            case 2
                plot(ax, rmsfe_tranq{m}, 'LineWidth',2, 'Color',colors(m,:), 'LineStyle','-');
            case 3
                plot(ax, rmsfe_covid{m}, 'LineWidth',2, 'Color',colors(m,:), 'LineStyle','-');
        end
    end
    xticks(ax, [1,3,5,7,9,11,13,15,17,19])
    xlim(ax, [0 20])
    if s==3
        ylim(ax, [0 10])
    else
        ylim(ax, [0 1.5])
    end
    set(ax, 'FontSize', 12)
    xticklabels(ax, {'135','120','110','95','85','75','60','50','35','15'})
    if s==1
        ylabel(ax, 'RMSFE','FontSize',16)
        title(ax, 'GFC','FontSize',16)
    elseif s==2
        title(ax, 'Tranquil','FontSize',16)
    else
        title(ax, 'Covid','FontSize',16)
    end
    hold(ax, 'off')
end

% Second row: CRPS
for s = 4:6
    ax = nexttile(s);
    hold(ax, 'on')
    for m = 1:nModels
        switch s
            case 4
                plot(ax, rtcrps_gfc{m}, 'LineWidth',2, 'Color',colors(m,:), 'LineStyle','-');
            case 5
                plot(ax, rtcrps_tranq{m}, 'LineWidth',2, 'Color',colors(m,:), 'LineStyle','-');
            case 6
                plot(ax, rtcrps_covid{m}, 'LineWidth',2, 'Color',colors(m,:), 'LineStyle','-');
        end
    end
    xticks(ax, [1,3,5,7,9,11,13,15,17,19])% TODO: double spacing of the available periods
    xlim(ax, [0 20])
    if s==6
        ylim(ax, [0 8])
    else
        ylim(ax, [0 1.5])
    end
    set(ax, 'FontSize', 12)
    xticklabels(ax, {'135','120','110','95','85','75','60','50','35','15'})
    if s==4
        ylabel(ax, 'CRPS','FontSize',16)
    end
    hold(ax, 'off')
end

% Add a legend spanning the top of the second row
lgd = legend(nexttile(2), model_names, 'NumColumns',2, 'Orientation','horizontal');
lgd.Layout.Tile = 'north'; % Place legend between rows

set(lgd,'FontSize',14);

% Add shared labels
xlabel(t, 'Days Until GDP Release','FontSize',16);

% Optional: adjust figure size/orientation as before
set(fig1,'PaperOrientation','landscape');
set(fig1,'PaperUnits','normalized');
set(fig1,'PaperPosition', [-0.05 0 1.1 1]);

if ~exist('Figures', 'dir')
    mkdir('Figures');
end

% Save as PDF and FIG
set(fig1, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);
t.TileSpacing = 'compact';
t.Padding = 'compact';
drawnow; % Ensure all graphics objects are rendered

exportgraphics(fig1, fullfile('Figures', 'nowcast_eval_comp.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'white', 'Resolution', 300);
saveas(fig1, fullfile('Figures', 'nowcast_eval_comp.fig'));

end

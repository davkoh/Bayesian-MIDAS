function trends = get_trends(mod, periods, y)
    % Extract trend data from the model
    for i = 1:length(periods)
        t = squeeze(mod.tau_all(:, periods(i), :));
      %  t = t - mean(t,2) +mean(y);
        trends(i).trend = t - mean(squeeze(median(t, 1))) + mean(y);
        trends(i).loc = median(trends(i).trend, 1); 
        trends(i).lower = quantile(trends(i).trend, 0.05); 
        trends(i).upper = quantile(trends(i).trend, 0.95);
    end
end
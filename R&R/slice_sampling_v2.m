% Function for slice sampling
function sample = slice_sampling_v2(current_value, log_likelihood,y_ar,x_ar, lower_bound, upper_bound)

    log_likelihood = @(x) -0.5*sum((y_ar - (2*x - 1).*x_ar).^2) + log(betapdf(x, 10, 2));

    % Set up initial interval
    y = log_likelihood(current_value) - exprnd(1);
    L = current_value - rand*abs(upper_bound - lower_bound);
    R = L + abs(upper_bound - lower_bound);
    
    % Iterate to find slice
    while log_likelihood(L) > y
        L = L - rand*abs(upper_bound - lower_bound);
    end
    while log_likelihood(R) > y
        R = R + rand*abs(upper_bound - lower_bound);
    end
    
    % Sample uniformly from slice
    while true
        sample = L + rand*(R - L);
        if log_likelihood(sample) > y
            break;
        elseif sample > current_value
            R = sample;
        elseif sample < current_value
            L = sample;
        else
            error('Slice sampler got stuck!');
        end
    end
end

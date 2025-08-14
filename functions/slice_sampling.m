% Function for slice sampling
function sample = slice_sampling(current_value, log_likelihood, lower_bound, upper_bound)
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

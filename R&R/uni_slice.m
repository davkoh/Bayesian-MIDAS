function x1 = uni_slice(x0, g, w, m, lower, upper, gx0)
    % Check the validity of the arguments
    if ~isnumeric(x0) || numel(x0)~=1 || ~isa(g, 'function_handle') || ...
            ~isnumeric(w) || numel(w)~=1 || w<=0 || ...
            ~isnumeric(m) || (~isinf(m) && (m<=0 || m>1e9 || floor(m)~=m)) || ...
            ~isnumeric(lower) || numel(lower)~=1 || x0<lower || ...
            ~isnumeric(upper) || numel(upper)~=1 || x0>upper || upper<=lower || ...
            (~isempty(gx0) && (~isnumeric(gx0) || numel(gx0)~=1))
        error('Invalid slice sampling argument');
    end

    % Keep track of the number of calls made to this function.
    % uni_slice_calls = uni_slice_calls + 1;

    % Find the log density at the initial point, if not already known.
    if isempty(gx0)
        gx0 = g(x0);
    end

    % Determine the slice level, in log terms.
    logy = gx0 - exprnd(1);

    % Find the initial interval to sample from.
    u = rand(1, 1)*w;
    L = x0 - u;
    R = x0 + (w - u);  % should guarantee that x0 is in [L,R], even with roundoff

    % Expand the interval until its ends are outside the slice, or until
    % the limit on steps is reached.
    if isinf(m)  % no limit on number of steps
        while true
            if L <= lower
                break;
            end
            if g(L) <= logy
                break;
            end
            L = L - w;
        end
        while true
            if R >= upper
                break;
            end
            if g(R) <= logy
                break;
            end
            R = R + w;
        end
    elseif m > 1  % limit on steps, bigger than one
        J = floor(rand(1, 1)*m);
        K = (m - 1) - J;
        while J > 0
            if L <= lower
                break;
            end
            if g(L) <= logy
                break;
            end
            L = L - w;
            J = J - 1;
        end
        while K > 0
            if R >= upper
                break;
            end
            if g(R) <= logy
                break;
            end
            R = R + w;
            K = K - 1;
        end
    end

    % Shrink interval to lower and upper bounds.
    if L < lower
        L = lower;
    end
    if R > upper
        R = upper;
    end

    % Sample from the interval, shrinking it on each rejection.
    while true
        x1 = L + rand*(R - L);
        gx1 = g(x1);
        if gx1 >= logy
            break;
        end
        if x1 > x0
            R = x1;
        else
            L = x1;
        end
    end

    
end

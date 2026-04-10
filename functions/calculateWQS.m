function wqs = calculateWQS(y_pred_all, yf, tau_grid, weight_type)
% CALCULATEWQS  Weighted Quantile Score from predictive draws
%   WQS = calculateWQS(Y_PRED_ALL, YF, TAU_GRID, WEIGHT_TYPE)
%
%   Inputs:
%     Y_PRED_ALL – MCMC × vint × nfor  array of predictive draws
%     YF         – nfor × 1  vector of actual outcomes
%     TAU_GRID   – vector of quantile levels (e.g. 0.05:0.05:0.95)
%     WEIGHT_TYPE – 1 = uniform, 2 = centre, 3 = tails
%
%   Output:
%     WQS – vint × nfor  matrix of weighted quantile scores

[~, vint, nfor] = size(y_pred_all);
ntau = length(tau_grid);

% Weight function
switch weight_type
    case 1, w = ones(1, ntau);                   % uniform
    case 2, w = tau_grid .* (1 - tau_grid);       % centre-weighted
    case 3, w = (2*tau_grid - 1).^2;              % tail-weighted
    otherwise, w = ones(1, ntau);
end
w = w / sum(w);

wqs = zeros(vint, nfor);
for t = 1:nfor
    yt = yf(t);
    for v = 1:vint
        draws = y_pred_all(:, v, t);
        q_hat = quantile(draws, tau_grid);          % 1 × ntau

        % Pinball (check) loss:  rho_tau(u) = u*(tau - I(u<0))
        u  = yt - q_hat;
        qs = u .* (tau_grid - (u < 0));             % non-negative scores

        wqs(v, t) = sum(w .* qs);
    end
end
end

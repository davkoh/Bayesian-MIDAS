function [dm_stat, pval] = dmtest_modified(e1, e2, h, loss_type)
% DMTEST_MODIFIED  Diebold-Mariano test with Harvey-Leybourne-Newbold correction
%   [DM, PV] = dmtest_modified(E1, E2) tests equal predictive accuracy
%   using squared loss.  DM > 0 means model 2 has smaller squared loss.
%   [DM, PV] = dmtest_modified(E1, E2, H) uses H-step-ahead HAC variance.
%   [DM, PV] = dmtest_modified(E1, E2, H, LOSS_TYPE) chooses the loss
%   differential. LOSS_TYPE = 'squared' uses E1.^2 - E2.^2,
%   LOSS_TYPE = 'raw' uses E1 - E2 directly.
%
%   Inputs:
%     E1, E2 – vectors of forecast errors (or loss values)
%     H      – forecast horizon for HAC variance (default: 1)
%     LOSS_TYPE – 'squared' (default) or 'raw'
%
%   Outputs:
%     DM_STAT – modified DM statistic (asympt. standard normal under H0)
%     PVAL    – two-sided p-value

if nargin < 3, h = 1; end
if nargin < 4, loss_type = 'squared'; end

e1 = e1(:);
e2 = e2(:);
T  = length(e1);

% Loss differential
if strcmp(loss_type, 'raw')
    d = e1 - e2;
else
    d = e1.^2 - e2.^2;
end
dbar = mean(d);

% Newey-West HAC variance with (h-1) lags
gamma0 = var(d, 1);            % population variance (1/T normalisation)
nw_sum = 0;
for j = 1:max(h-1, 0)
    gj     = (d(j+1:end) - dbar)' * (d(1:end-j) - dbar) / T;
    nw_sum = nw_sum + 2 * gj;
end
V = (gamma0 + nw_sum) / T;

if V <= 0
    dm_stat = 0;
    pval    = 1;
    return
end

dm_raw = dbar / sqrt(V);

% HLN small-sample correction (Harvey, Leybourne & Newbold, 1997)
hlncorr = sqrt((T + 1 - 2*h + h*(h-1)/T) / T);
dm_stat = dm_raw * hlncorr;

pval = 2 * (1 - normcdf(abs(dm_stat)));
end

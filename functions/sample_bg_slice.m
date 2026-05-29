function [bg] = sample_bg_slice(bg, varphis, bg_lb, bg_ub, c, d)
%SAMPLE_BG_SLICE Slice-sample b_g | {varphi_{g,i}^2} under
%   varphi_{g,i}^2 | b_g ~ Gamma(b_g, 1) (shape, rate),
%   b_g           ~ Gamma(c, d)        (shape, rate).
%
%   bg      : current value of b_g
%   varphis : vector of within-group local variances for group g
%   bg_lb   : lower bound for slice support
%   bg_ub   : upper bound for slice support
%   c, d    : hyperparameters of the Gamma prior on b_g

n              = numel(varphis);
sum_log_varphi = sum(log(varphis));
sum_varphi     = sum(varphis);

log_likelihood = @(x) -n*log(gamma(x)) + (x-1)*sum_log_varphi - sum_varphi ...
                      + c*log(d) - log(gamma(c)) + (c-1)*log(x) - d*x;

bg = uni_slice(bg, log_likelihood, 1, inf, bg_lb, bg_ub, []);

end

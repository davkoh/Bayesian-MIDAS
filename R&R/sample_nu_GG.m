function [nu f_nu] = sample_nu_GG(lam,nu_ub,n_grid)
T = size(lam,1);
sum1 = sum(log(lam));
sum2 = sum(1./lam);
f_nu = @(x) T*(x/2.*log(x/2)-gammaln(x/2)) ...
- (x/2+1)*sum1 - x/2*sum2;
nu_grid = linspace(2+rand/100,nu_ub-rand/100,n_grid)’;
lp_nu = f_nu(nu_grid); % log-density of nu
p_nu = exp(lp_nu - max(lp_nu)); % density of nu (unnormalized)
p_nu = p_nu/sum(p_nu); % density of nu (normalized)
cdf_nu = cumsum(p_nu); % cdf of nu
nu = nu_grid(find(rand<cdf_nu, 1));
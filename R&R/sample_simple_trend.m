function [output_trend] = sample_simple_trend(ystar,prior_trend)

T = length(ystar);

% Unpack data
sig2_y = prior_trend.sig2_y;
sig2_tau = prior_trend.sig2_tau;
tau0 = prior_trend.tau0;
H = prior_trend.H;
nu_omega0 = prior_trend.nu_omega0;
S_omega0 = prior_trend.S_omega0;
a0 = prior_trend.a0;
b0 = prior_trend.b0;
tau = prior_trend.tau;

% sample tau
Ktau = H'*H/sig2_tau + speye(T)./sig2_y;
Ctau = chol(Ktau,"lower");
tau_hat = Ktau\(tau0/sig2_tau*H'*H*ones(T,1)+ speye(T)./sig2_y*ystar);
tau = tau_hat + Ctau'\randn(T,1);

% Sample omega
sig2_tau = 1/gamrnd(nu_omega0 + T/2, ...
1/(S_omega0 + (tau-tau0)'*H'*H*(tau-tau0)/2));


% sample tau0
Ktau0 = 1/b0 + 1/sig2_tau;
tau0_hat = Ktau0\(a0/b0 + tau(1)/sig2_tau);
tau0 = tau0_hat + sqrt(Ktau0)'\randn;

output_trend.sig2_y = sig2_y;
output_trend.sig2_tau = sig2_tau;
output_trend.tau0 = tau0;
output_trend.H = H;
output_trend.nu_omega0 = nu_omega0;
output_trend.S_omega0 = S_omega0;
output_trend.a0 = a0;
output_trend.b0 = b0;
output_trend.tau = tau;



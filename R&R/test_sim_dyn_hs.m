%% The Local Trend Model with SV with Horseshoe Prior %%

clear all;

% Choose Parameters of the DGP
T = 200;
sd_tau = 0.01;
t_0 = 0;
sd_y = 0.1;

% Generate Data
tau = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);

t = [0.1:0.1:20]
a = sin(t);

for t= 1:T
    if t == 1
        tau(1) = t_0 + sd_tau*randn;
    else
        tau(t) = a(t-1) + sd_tau*randn;
    end

    if t == 100
        tau(t) = a(t-1) + sd_tau*randn;
        y(t) = tau(t) + 10 + sd_y*randn;
    else
        y(t) = tau(t) + sd_y*randn;
    end
end


tau_true = tau;
h_true = h;

clf;
plot(y)
plot(tau_true)

% Priors
nu_sig0 = 0.1; S_sig0 = 0.1;


% Initalise series
tau0 = 0;
lam2t = zeros(T,1);
nu2 = 0;
phi = 0.9;
Hphi = speye(T) - sparse(2:T,1:(T-1),phi*ones(1,T-1),T,T);
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
mu = 0;
mu_tilde = [mu;(1-phi)*ones(T-1,1)*mu];
h_tilde = h-mu;
xi = ones(T,1); % diagonals of the prior on h_tilde
iSig_xi = sparse(1:T,1:T,1./xi);
v = ones(T,1); % diagonals of mixture rep of the likelihood
iSig_v = sparse(1:T,1:T,1./v);
h = zeros(T,1);
xi_0 = 0;
xi_mu = 0;
sig2_y = 1;

% Sampler Specs
rng(1,'twister');
iter = 10000;
burn = iter/2;
mcmc = iter-burn;

% Function for log posterior of phi
log_likelihood = @(x) -0.5*sum((y_ar - (2*x - 1).*x_ar).^2) + log(betapdf(x, prior_phi(1), prior_phi(2)));    
prior_phi = [1, 1]; % Parameters for Beta prior

% Prep input the sampling function
    input_dyn_hs_trend.tau = tau;
    input_dyn_hs_trend.h = h;
    input_dyn_hs_trend.phi = phi;
    input_dyn_hs_trend.h_tilde = h_tilde;
    input_dyn_hs_trend.Hphi = Hphi;
    input_dyn_hs_trend.iSig_xi = iSig_xi; 
    input_dyn_hs_trend.mu_tilde = mu_tilde;
    input_dyn_hs_trend.mu = mu;
    input_dyn_hs_trend.xi = xi;
    input_dyn_hs_trend.xi_mu = xi_mu;
    input_dyn_hs_trend.xi_0 = xi_0;
    input_dyn_hs_trend.sig2_y = sig2_y;
    input_dyn_hs_trend.y = y;
    input_dyn_hs_trend.prior_phi = prior_phi;
    input_dyn_hs_trend.H = H;
    
%% Gibbs Sampler

% Storage
tau_save = NaN(T,mcmc);
phi_save = NaN(mcmc,1);
sig2_save = NaN(mcmc,1);
h_save = NaN(T,mcmc);
mu_save = NaN(mcmc,1);

tic
for i = 1:iter
    %% 1: Sample Trend
        ystar = log([tau(1);tau(2:end)-tau(1:end-1)].^2 + 0.0001);
        [input_dyn_hs_trend] = sample_dyn_hs_trend_v2(ystar,input_dyn_hs_trend);

        tau = input_dyn_hs_trend.tau;

    %% 2: Sample observation variance
        sig2_y = 1/gamrnd(nu_sig0 + T/2,1/(S_sig0 + (y-tau)'*(y-tau)/2));
       input_dyn_hs_trend.sig2_y = sig2_y;

if i> burn 
    isave = i - burn;
    tau_save(:,isave) = tau;
    phi_save(isave) = input_dyn_hs_trend.phi;
    sig2_save(isave) = sig2_y;
    h_save(:,isave) = input_dyn_hs_trend.h;
    mu_save(isave) = input_dyn_hs_trend.mu;

    
end

end
toc

ptau = mean(tau_save,2);
clf;
plot(ptau)
hold on
plot(y)
plot(tau_true)
plot(exp(mean(h_save,2)))

sqrt(mean(sig2_save))
sqrt(sum((ptau-y).^2)/T)
sqrt(sum(y-tau_true).^2/T)
sqrt(sum(ptau(1:50)-tau_true(1:50)).^2/50)


ph = mean(h_save,2);
clf;
plot(exp(0.5*ph))
hold on
plot(exp(0.5*h_true))

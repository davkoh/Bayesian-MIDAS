%% The Local Trend Model with SV with Horseshoe Prior in both + t-errors %%

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

t = [0.1:0.1:T/10]
a = sin(t);

for t= 1:T
    if t == 1
        tau(1) = t_0 + sd_tau*randn;
    else
        tau(t) = a(t-1) + sd_tau*randn;
    end

    if t == 100
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
nu_sig0 = 3; S_sig0 = 1*(nu_sig0-1);


% Initalise series
tau0 = 0;
lam2t = zeros(T,1);
nu2 = 0;
phi = 0; % Start at zero perhaps
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
lamt = ones(T,1);
nut = 5;
alpha = 2;
beta = .2;
nu_lb = 2;
nu_ub = 30; 
sig2_tau = .1; 
tau0 = 0;
a0 = 5; b0 = 100;
nu_sig0 = 3; S_sig0 = 1*(nu_sig0-1);
nu_omega0 = 3; S_omega0 = .25^2*(nu_omega0-1);

% Sampler Specs
rng(1,'twister');
iter = 10000;
burn = iter/2;
mcmc = iter-burn;

% Function for log posterior of phi
log_likelihood = @(x) -0.5*sum((y_ar - (2*x - 1).*x_ar).^2) + log(betapdf(x, prior_phi(1), prior_phi(2)));    
prior_phi = [10, 2]; % Parameters for Beta prior

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

    input_dyn_hs_sv.g = h;
    input_dyn_hs_sv.phi = phi;
    input_dyn_hs_sv.g_tilde = h_tilde;
    input_dyn_hs_sv.Hphi = Hphi;
    input_dyn_hs_sv.iSig_xi = iSig_xi; 
    input_dyn_hs_sv.mu_tilde = mu_tilde;
    input_dyn_hs_sv.mu = mu;
    input_dyn_hs_sv.xi_mu = xi_mu;
    input_dyn_hs_sv.xi_0 = xi_0;
    input_dyn_hs_sv.y = y;
    input_dyn_hs_sv.prior_phi = prior_phi; % Parameters for Beta prior
    input_dyn_hs_sv.H = H;
    

    input_trend.sig2_y = sig2_y;
    input_trend.sig2_tau = sig2_tau;
    input_trend.tau0 = tau0;
    input_trend.H = H;
    input_trend.nu_omega0 = nu_omega0;
    input_trend.S_omega0 = S_omega0;
    input_trend.a0 = a0;
    input_trend.b0 = b0;
    input_trend.tau = tau;
%% Gibbs Sampler

% Storage
tau_save = NaN(T,mcmc);
phih_save = NaN(mcmc,1);
phitau_save = NaN(mcmc,1);
h_save = NaN(T,mcmc);
muh_save = NaN(mcmc,1);
mutau_save = NaN(mcmc,1);
nu_save = NaN(mcmc,1);

tic
for i = 1:iter
    %% 1: Sample Trend
        ystar = y;
        [input_trend] = sample_simple_trend(ystar,input_trend);

        tau = input_trend.tau;

    %% 2: Sample observation variance
        ystar = log((y-tau).^2./lamt + 0.0001);
        [input_dyn_hs_sv] = sample_dyn_hs_sv(ystar,input_dyn_hs_sv);
        input_trend.sig2_y = exp(input_dyn_hs_sv.g).*lamt ;

    %% 3.1: Sample Mixture
    lamt = 1./gamrnd((nut+1)/2,2./(nut+(1./exp(input_dyn_hs_sv.g).*speye(T))*((y-tau).^2)));

    %% 3.2: Sample degrees of freedom
    nut = sample_nu_slice(nut,lamt,nu_lb,nu_ub,alpha,beta);
    input_dyn_hs_trend.sig2_y = exp(input_dyn_hs_sv.g).*lamt ;

if i> burn 
    isave = i - burn;
    tau_save(:,isave) = tau;
    phitau_save(isave) = input_dyn_hs_trend.phi;
    phih_save(isave) = input_dyn_hs_sv.phi;
    h_save(:,isave) = lamt;%input_dyn_hs_sv.g;
    muh_save(isave) = input_dyn_hs_sv.mu;
    mutau_save(isave) = input_dyn_hs_trend.mu;
    nu_save(isave) = nut;
end

end
toc

ptau = mean(tau_save,2);
clf;
plot(ptau)
hold on
plot(tau_true)
plot(y)
ph = mean(h_save,2);
plot(ph)

sqrt(sum((ptau-y).^2)/T)
sqrt(sum((ptau-tau_true).^2)/T)

ph = mean(h_save,2);
clf;
plot(exp(1/2*ph))
hold on
plot(exp(0.5*h_true))

clf;
histogram(nu_save)
plot(exp(h_save(140,:)))

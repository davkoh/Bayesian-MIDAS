%% The Local Trend Model with SV with Horseshoe Prior %%

clear all;

% Choose Parameters of the DGP
T = 70;
sd_tau = 0.1;
t_0 = 0;
sd_y = 1;
h_0 = 0;
sd_h = 1;

% Generate Data
tau = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);

outlier = 150;


for t= 1:T
    if t == 1
        tau(1) = t_0 + sd_tau*randn;
        h(1) = h_0 + sd_h*randn;
    else
    tau(t) = tau(t-1) + sd_tau*randn;
    h(t) = h(t-1) + sd_h*randn;
    end

    
        y(t) = tau(t) + exp(1/2*h(t))*trnd(2);
    
end


tau_true = tau;
h_true = h;

clf;
plot(y)
plot(h_true)
plot(tau_true)

% Priors
V_t0 = 10; % Variance of initial condition
Vh = 10;

prior_tau.V_t0 = 10;
prior_tau.H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);


prior_sv.Vh = 10;



% Initalise series
lambdatau = ones(T,1);
nutau = 1;
eta_lam = ones(T,1);
eta_nu = 1;
tau0 = 0;
Sigmatau = ones(T,1);
sig2 = 1;
nuh=1;
lambdah=ones(T,1);
Sigmah=(lambdah.*nuh);
h = ones(T,1);
h0 = log(var(y));
tau= zeros(T,1);
nub = 50;
lamt = ones(T,1);
nut = 5;

% Sampler Specs
rng(1,'twister');
iter = 20000;
burn = iter/2;
mcmc = iter-burn;

%% Gibbs Sampler

% Storage
tau_save = NaN(T,mcmc);
tau0_save = NaN(mcmc,1);
sig2_save = NaN(mcmc,1);
lam_save = NaN(T,mcmc);
nu_save = NaN(mcmc,1);
h_save = NaN(T,mcmc);
h0_save = NaN(mcmc,1);
nut_save = NaN(mcmc,1);

for i = 1:iter
    %% Sample tau
    [tau,Sigmatau,tau0,lambdatau,nutau] = sample_trend_unitroot_hs(y,tau,Sigmatau,tau0,lambdatau,nutau,exp(h).*lamt,prior_tau);

    %% Sample h
    u = y-tau;
    Ystar = log(u.^2.*lamt.^(-1) + .0001);
    [h,Sigmah,h0,lambdah,nuh] = sample_sv_unitroot_hs(Ystar,h,Sigmah,h0,lambdah,nuh,prior_sv);

    %% Sample t-errors
    [lamt,nut] = sample_t_errors(u,lamt,nut,nub,exp(h));
    

if i> burn 
    isave = i - burn;
    tau_save(:,isave) = tau;
    tau0_save(isave) = tau0;
    lam_save(:,isave) = lambdatau;
    nu_save(isave) = nutau;
    h_save(:,isave) = h;
    h0_save(isave) = h0;
    nut_save(isave) = nut;
end

end

ptau = mean(tau_save,2);
clf;
plot(ptau)
hold on
plot(tau_true)

ph = mean(h_save,2);
clf;
plot(exp(0.5*ph))
hold on
plot(exp(0.5*h_true))

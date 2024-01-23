%% The Local Trend Model with SV, Horseshoe Prior  in SV%%

clear all;

% Choose Parameters of the DGP
T = 80;
sd_tau = 0.1;
t_0 = 0;
sd_y = 1;
h_0 = 0;
sd_h = 0.1;

% Generate Data
tau = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);


for t= 1:T
    if t == 1
        tau(1) = t_0 + sd_tau*randn;
        h(1) = h_0 + sd_h*randn;
    else
    tau(t) = tau(t-1) + sd_tau*randn;
    h(t) = h(t-1) + sd_h*randn;
    end

    if t == 70
        tau(t) = tau(t);
        h(t) = h(t);
        y(t) = tau(t) + exp(1/2*h(t))*randn -20;
    elseif t ==71
        tau(t) = tau(t);
        h(t) = h(t);
        y(t) = tau(t) + exp(1/2*h(t))*randn +20;
    else
        y(t) = tau(t) + exp(1/2*h(t))*randn;
    end
end


tau_true = tau;
h_true = h;

clf;
plot(y)
plot(h_true)
plot(tau_true)

% Priors
V_t0 = 10; % Variance of initial condition
nu_sig0 = 3; S_sig0 = 1*(nu_sig0-1); % \sigma^2
Vh = 10;
nu_omega0 = 3; S_omega0 = .25^2*(nu_omega0-1); % \omega_tau^2

% Set up difference matrix
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);

% Initalise series
lam = ones(T,1);
nu = 1;
eta_lam = ones(T,1);
eta_nu = 1;
tau0 = 0;
Sigma_tau = ones(T,1);
sig2 = 1;
tauh=1;
lambdah=ones(T,1);
Sigh=(lambdah.*tauh);
h = ones(T,1);
h0 = log(var(y));

 omega_tau2 = .1;

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

for i = 1:iter
    % Sample \tau
    K_tau = H'*H/omega_tau2 +  1./exp(h).*speye(T);
    C_tau = chol(K_tau,"lower");
    tau_hat = K_tau\(tau0*H'*H*ones(T,1)/omega_tau2 + 1./exp(h).*speye(T)*y);
    tau = tau_hat + C_tau'\randn(T,1);

    % sample omega2
    omega_tau2 = 1/gamrnd(nu_omega0 + T/2, ...
    1/(S_omega0 + (tau-tau0)'*H'*H*(tau-tau0)/2));

    % Sample \tau_0
    Ktau0 = 1/V_t0 + 1/omega_tau2;
    tau0_hat = Ktau0\(tau(1)/omega_tau2);
    tau0 = tau0_hat + sqrt(Ktau0)'\randn;

    %% Sample h
    u = y-tau;
    Ystar = log(u.^2 + .0001);
    h = SVRW2(Ystar,h,Sigh,h0);

    %% sample h0
    Kh0 = 1./Sigh(1)+ 1./Vh;
    h0hat = Kh0\(h(1)./Sigh(1));
    h0 = h0hat + chol(Kh0,'lower')'\randn;

     %% sample Sigh:
    
        e = (h - [h0;h(1:T-1,:)]);
        v=1./gamrnd(1,1./(1+1./lambdah));
        lambdah=1./gamrnd(1, 1./(1./v + (e.^2)./(2.*tauh)));
        xi=1./gamrnd(1,1./(1+1./tauh));
        tauh=1./gamrnd((T+1)/2, 1./(1./xi +0.5*sum(sum(e.^2./lambdah))  ));
        Sigh=(lambdah.*tauh)+1e-10;
    

if i> burn 
    isave = i - burn;
    tau_save(:,isave) = tau;
    tau0_save(isave) = tau0;
    sig2_save(isave) = sig2;
    lam_save(:,isave) = lam;
    nu_save(isave) = nu;
    h_save(:,isave) = h;
    h0_save(isave) = h0;
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

plam = mean(lam_save,2);

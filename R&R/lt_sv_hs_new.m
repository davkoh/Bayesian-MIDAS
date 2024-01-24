%% T-SV with new sampler %%
clear all;

% Choose Parameters of the DGP
T = 200;
sd_tau = 0.1;
t_0 = 0;
sd_y = 1;
h_0 = 2;
sd_h = .5;
phih = 0.98;

% Generate Data
tau = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);
tau = zeros(T,1);
outlier = 20;


for t= 1:T
    if t == 1
        tau(1) = t_0 + sd_tau*randn;
        h(1) = h_0 + sd_h*randn;
    else
    tau(t) = tau(t-1) + sd_tau*randn;
    h(t) = phih*h(t-1) + sd_h*randn;
    end

    if t == 70
        tau(t) = tau(t);
        h(t) = h(t);
        y(t) = tau(t) + exp(1/2*h(t))*randn -outlier;
    elseif t ==71
        tau(t) = tau(t);
        h(t) = h(t);
        y(t) = tau(t) + exp(1/2*h(t))*randn +outlier;
    else
        y(t) = tau(t) + exp(1/2*h(t))*randn;
    end
end



h_true = h;
tau_true = tau;

clf;
plot(y)
plot(h_true)
plot(tau)


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
omegah2 = 0.1;
muh = 0;
lamh = ones(T,1);
tauh = 1;
Sigmah = ones(T,1);

 omega_tau2 = .1;

% Sampler Specs
rng(1,'twister');
iter = 20000;
burn = iter-iter/2;
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

mu_save = NaN(mcmc,1);
phih_save = NaN(mcmc,1);
omegah_save = NaN(mcmc,1);



sv_prior.pr_S_sigh = 0.1; % Prior on the state variation;
sv_prior.pr_nu_sigh = 5;

sv_prior.phih0 = .95; sv_prior.Vphih = .5^2;

sv_prior.muh0 = 1; sv_prior.Vmuh = 10;




for i = 1:iter

     % Sample \tau
    K_tau = H'*sparse(diag(1./Sigma_tau))*H +  1./exp(h).*speye(T);
    C_tau = chol(K_tau,"lower");
    tau_hat = K_tau\(tau0*H'*sparse(diag(1./Sigma_tau))*H*ones(T,1) + 1./exp(h).*speye(T)*y);
    tau = tau_hat + C_tau'\randn(T,1);

    % Sample \Sigma_tau
        e = tau-[tau0;tau(1:T-1)];
        lam = 1./gamrnd(1, 1./( 1./eta_lam + 0.5*e.^2/nu));
        % sample Global
        nu = 1/gamrnd( 0.5*T, 1/( 1/eta_nu + 0.5*sum(sum(e.^2./lam))  ) );
        % sample mixture local
        eta_lam = 1./gamrnd(1, 1./(1 + 1./lam));
        % samplel mixture
        eta_nu = 1/gamrnd(1, 1/( 1 + 1/nu));
    Sigma_tau = nu*lam;

    % Sample \tau_0
    Ktau0 = 1/V_t0 + 1/(Sigma_tau(1));
    tau0_hat = Ktau0\(tau(1)/Sigma_tau(1));
    tau0 = tau0_hat + sqrt(Ktau0)'\randn;

    % Sample h

s2 = (y-tau).^2;

[h,phih,Sigmah,muh,lamh,tauh] = sample_sv_hs_ar1(s2,h,sv_prior,phih,Sigmah,muh,lamh,tauh);

    

    if i> burn 
    isave = i - burn;
    omegah_save(isave) = sqrt(omegah2);
    h_save(:,isave) = h;
    mu_save(isave) = muh;
    phih_save(isave) = phih;
    tau_save(:,isave) = tau;
    tau0_save(isave) = tau0;
end


end



% Evaluation

ptau = mean(tau_save(:,1:1:end),2);
clf;
plot(ptau)
hold on
plot(tau_true)

ph = mean(h_save(:,1:1:end),2);
clf;
plot(exp(0.5*ph))
hold on
plot(exp(0.5*h_true))


hhat = mean(exp(h_save/2),2);  %% plot std dev
plot(hhat)

phihat = phih_save(1:100:end,:);



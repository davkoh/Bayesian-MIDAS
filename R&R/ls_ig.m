%% The Local Trend Model with IG Priors %%

% Choose Parameters of the DGP
T = 200;
sd_tau = .5;
t_0 = 0;
sd_y = 0.5;

% Generate Data
tau = zeros(T,1);
y = zeros(T,1);

for t= 1:T
    if t == 1
        tau(1) = t_0 + sd_tau*randn;
    else
    tau(t) = tau(t-1) + sd_tau*randn;
    end

    if t == 100
        tau(t) = tau(t) + 10;
        y(t) = tau(t) + sd_y*randn;
    else
        y(t) = tau(t) + sd_y*randn;
    end
end

tau_true = tau;

% Priors
V_t0 = 10; % Variance of initial condition
nu_sig0 = 3; S_sig0 = 1*(nu_sig0-1); % \sigma^2

% Set up difference matrix
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);

% Initalise series
lam = ones(T,1);
nu = 1;
eta_lam = ones(T,1);
eta_nu = 1;
tau0 = 0;
Sigma_tau = 1;
sig2 = 1;

% Sampler Specs
rng(1,'twister');
iter = 8000;
burn = iter/2;
mcmc = iter-burn;

%% Gibbs Sampler

% Storage
tau_save = NaN(T,mcmc);
tau0_save = NaN(mcmc,1);
sig2_save = NaN(mcmc,1);
lam_save = NaN(T,mcmc);
nu_save = NaN(mcmc,1);

for i = 1:iter
    % Sample \tau
    K_tau = H'*H/Sigma_tau + 1/sig2*speye(T);
    C_tau = chol(K_tau,"lower");
    tau_hat = K_tau\(tau0/Sigma_tau*(H)'*H*ones(T,1) + y/sig2);
    tau = tau_hat + C_tau'\randn(T,1);
    

    Sigma_tau = 1/gamrnd(5 + T/2, 1/(.1^2*(5-1) + (tau-tau0)'*(H)'*H*(tau-tau0)/2));

    % Sample \tau_0
    Ktau0 = 1/V_t0 + 1/(Sigma_tau(1));
    tau0_hat = Ktau0\(tau(1)/Sigma_tau(1));
    tau0 = tau0_hat + sqrt(Ktau0)'\randn;

    % Sample sig2
    sig2 = 1/gamrnd(5 + T/2,1/(S_sig0 + (y-tau)'*(y-tau)/2));

if i> burn 
    isave = i - burn;
    tau_save(:,isave) = tau;
    tau0_save(isave) = tau0;
    sig2_save(isave) = sig2;
    %lam_save(:,isave) = lam;
    %nu_save(isave) = nu;
end

end

ptau = mean(tau_save,2);
clf;
plot(ptau,LineWidth=2)
hold on
plot(tau_true)

plam = mean(lam_save,2);

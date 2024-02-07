%% Testing the Non-centred hierarchical horseshoe prior Model %%

clear all;
rng(1,'twister');


% Generate Data

% Choose Parameters of the DGP
T = 200;
sd_y = 0.1;
alpha = 0;
omega_tau = 0.2;



% Generate Data
tau = zeros(T,1);
tautilde = zeros(T,1);
y = zeros(T,1);

for t= 2:T
tautilde(t) = tautilde(t-1) + randn;  
y(t) = alpha + omega_tau * tautilde(t) + sd_y*randn;

end

tau_true = alpha + omega_tau*tautilde;
tautilde_true = tautilde;
omegatau_true = omega_tau;

clf;
plot(y)
plot(tau_true)
hold on
plot(tautilde_true)


%% Priors
% sigma_y^2 ~ IG(a,b)
nu_sig0 = 3; S_sig0 = 1*(nu_sig0-1);
% \tilde{tau}_0 ~ N(0,V_tau^2)
nu_tautilde0 = 0.01; S_tautilde0 = 0.01;

% Pre-compute useful things
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Hinv = H\speye(T);
Pgam = H'*H;

%% Starting Values
beta =  0.1*ones(2,1);
lamba_beta = ones(2,1);
nu_beta = 1;
eta_lam = ones(2,1);
eta_nubeta = 1;
tautilde_0 = 0;
sigy2 = 0.1;
V_tautilde0 = 0.1;
Sigma_tau = nu_beta*lamba_beta;

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
    %% Sample tautilde
    Xgam = omega_tau*speye(T);
    Kgam = Pgam + Xgam'*Xgam/sigy2;
    gam_hat = Kgam\( Hinv*tautilde_0*ones(T,1) + 1/sigy2*Xgam'*(y-alpha));
    tau = gam_hat + chol(Kgam,'lower')'\randn(T,1);
   

    %% Sample beta (alpha,omegatau)
    X = [ones(T,1) tau];
    Kbeta = diag(1./Sigma_tau) + 1/sigy2 * X'*X;
    beta_hat = Kbeta\(1/sigy2*X'*y);
    beta = beta_hat + chol(Kbeta,'lower')'\randn(2,1);
    alpha = beta(1); omegatau = beta(2);

    %% permutate the signs of tau and sigtau
    U = -1 + 2*(rand>0.5);
    tau = U*tau;
    beta(2) = U*beta(2);
    alpha = beta(1); omegatau = beta(2);
    tautilde_0 = U*tautilde_0;

    %% Sample Sigma_tau
        % Sample 
        lam_beta = 1./gamrnd(1, 1./( 1./eta_lam + 0.5*beta.^2/nu_beta));
        % sample Global
        nu_beta = 1/gamrnd( 0.5*2, 1/( 1/eta_nubeta + 0.5*sum(sum(beta.^2./lam_beta))  ) );
        % sample mixture local
        eta_lam = 1./gamrnd(1, 1./(1 + 1./lam_beta));
        % samplel mixture
        eta_nu = 1/gamrnd(1, 1/( 1 + 1/nu_beta));
    Sigma_tau = nu_beta*lam_beta;

    %% Sample tautilde_0
    Ktautilde_0 = 1/(V_tautilde0) + 1;
    tautilde_0_hat = Ktautilde_0\(tau(1)/(V_tautilde0));
    tautilde_0 = tautilde_0_hat + chol(Ktautilde_0,'lower')'\randn;

    %% Sample V_tautilde0
    V_tautilde0 = 1/gamrnd(nu_tautilde0 + 1/2,1/(S_tautilde0 + (tau(1)-alpha)'*(tau(1)-alpha)/2));

    %% Sample sigy2
    sigy2 = 1/gamrnd(nu_sig0 + T/2,1/(S_sig0 + (y-alpha - omegatau*tau)'*(y-alpha - omegatau*tau)/2));






end


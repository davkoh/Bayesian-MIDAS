function [tau,Sigmatau,tau0,lambdatau,nutau] = sample_trend_unitroot_hs(y,tau,Sigmatau,tau0,lambdatau,nutau,diag_obs,prior_tau)

%% Function Description

 % Unpack prior
 V_t0 = prior_tau.V_t0;
 H = prior_tau.H;
 T = length(diag_obs);

 % Sample \tau
K_tau = H'*sparse(diag(1./Sigmatau))*H +  1./diag_obs.*speye(T);
C_tau = chol(K_tau,"lower");
tau_hat = K_tau\(tau0*H'*sparse(diag(1./Sigmatau))*H*ones(T,1) + 1./diag_obs.*speye(T)*y);
tau = tau_hat + C_tau'\randn(T,1);

% Sample \Sigmatau
    e = tau-[tau0;tau(1:T-1)];
    % sample mixture local
    eta_lam = 1./gamrnd(1, 1./(1 + 1./lambdatau));
    lambdatau = 1./gamrnd(1, 1./( 1./eta_lam + 0.5*e.^2/nutau));
    % samplel mixture
    eta_nu = 1/gamrnd(1, 1/( 1 + 1/nutau));
    % sample Global
    nutau = 1/gamrnd( 0.5*T, 1/( 1/eta_nu + 0.5*sum(sum(e.^2./lambdatau))  ) );

Sigmatau = nutau*lambdatau;

% Sample \tau_0
Ktau0 = 1/V_t0 + 1/(Sigmatau(1));
tau0_hat = Ktau0\(tau(1)/Sigmatau(1));
tau0 = tau0_hat + sqrt(Ktau0)'\randn;
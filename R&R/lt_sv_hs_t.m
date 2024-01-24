%% The Local Trend Model with SV with Horseshoe Prior %%

clear all;

% Choose Parameters of the DGP
T = 200;
sd_tau = 0.1;
t_0 = 0;
sd_y = 1;
h_0 = 0;
sd_h = 0.1;

% Generate Data
tau = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);

outlier = 20;


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
        y(t) = tau(t) + exp(1/2*h(t))*randn -outlier;
    elseif t ==71
        tau(t) = tau(t);
        h(t) = h(t);
        y(t) = tau(t) + exp(1/2*h(t))*randn +outlier;
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
nut = 5;
lamt= ones(T,1);

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
lamt_save = NaN(T,mcmc);
nut_save = NaN(mcmc,1);

for i = 1:iter
    % Sample \tau
    K_tau = H'*sparse(diag(1./Sigma_tau))*H +  1./exp(h).*1./lamt.*speye(T);
    C_tau = chol(K_tau,"lower");
    tau_hat = K_tau\(tau0*H'*sparse(diag(1./Sigma_tau))*H*ones(T,1) + 1./exp(h).*1./lamt.*speye(T)*y);
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

    %% Sample h
    u = (y-tau);
    Ystar = log((u.^2).*lamt.^(-1) + .0001);
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

        %% Sample t-errors
       % h = zeros(T,1);
lamt = 1./gamrnd((nut+1)/2,2./(nut+(1./exp(h).*speye(T))*(u.^2)));

% sample nu
[nut,flag] = sample_nu_MH(lamt,nut,50);


if i> burn 
    isave = i - burn;
    tau_save(:,isave) = tau;
    tau0_save(isave) = tau0;
    sig2_save(isave) = sig2;
    lam_save(:,isave) = lam;
    nu_save(isave) = nu;
    h_save(:,isave) = h;
    h0_save(isave) = h0;
    lamt_save(:,isave) = lamt;
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

plam = mean(lam_save,2);

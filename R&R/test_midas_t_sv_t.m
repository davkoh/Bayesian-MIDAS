%% The Local Trend Model with SV with Horseshoe Prior %%

clear all;

% Choose Parameters of the DGP
    % Time-series
    T = 100;
    sd_tau = 0.5;
    t_0 = 0;
    sd_y = 1;
    h_0 = 0;
    sd_h = .2;
    % Group-correlated Data
    G = 5;
    K = 50;
    corr_x = 0.8;
    Sigma_X = kron(eye(G),ones(K/G)*corr_x);
    Sigma_X = Sigma_X- diag(diag(Sigma_X)) + diag(ones(K,1));
    X = mvnrnd(zeros(K,1),Sigma_X,T);
    beta = [G/G*ones(G,1);zeros(G,1);zeros(G,1);zeros(G,1);G/G*(ones(G,1));zeros(G,1);zeros(G,1);zeros(G,1);zeros(G,1);zeros(G,1)];
    beta_true = beta;


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

        if (t == 70)
        tau(t) = tau(t);
        h(t) = h(t);
        y(t) = X(t,:)*beta + tau(t) + exp(1/2*h(t))*randn -outlier;
    elseif (t ==71)
        tau(t) = tau(t);
        h(t) = h(t);
        y(t) = X(t,:)*beta + tau(t) + exp(1/2*h(t))*randn +outlier;
    else
        y(t) = X(t,:)*beta + tau(t) + exp(1/2*h(t))*randn;
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

% Set up data for grouped prior
    grp_midas = [];

    for i = 1:G
        grp_midas = [grp_midas;ones(K/G,1)*i];
    end
    
hyperpars = [1/T,0.5];
sum_grp = size(unique(grp_midas),1);

% Data
data.Y = y;
data.X = X;
data.tX = X';
data.gl_param_expand_diag_inv = zeros(K,1);
data.K = K;
data.tau_shape_const = (K+1)/2;
data.stable_const = 1e-07;
data.G = G;
data.grp_idx = grp_midas;
data.grp_size = histc(grp_midas, unique(grp_midas));
data.grp_size_cs = cumsum(histc(grp_midas, unique(grp_midas)));
data.a = repmat(1/T,sum_grp,1);
data.b = repmat(0.5,sum_grp,1);
data.iOh = eye(T);

% Hyperparameters
hyper_params.tau_sq = 1;
hyper_params.nu = 5;
hyper_params.gamma_sq = ones(G,1);
hyper_params.lambda_sq = ones(K,1);


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
beta_save = NaN(K,mcmc);

for i = 1:iter

    %% Sample beta
    data.Y = y-tau;
    data.iOh = 1./(exp(h).*lamt).*speye(T);
    [beta,hyper_params,data] = sample_beta(data,hyper_params);

    %% Sample tau
    [tau,Sigmatau,tau0,lambdatau,nutau] = sample_trend_unitroot_hs(y-X*beta,tau,Sigmatau,tau0,lambdatau,nutau,exp(h).*lamt,prior_tau);

    %% Sample h
    u = y-tau-X*beta;
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
    beta_save(:,isave) = beta;
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

pbeta = mean(beta_save,2);
clf;
plot(pbeta)
hold on
plot(beta_true)


clf;
hist(nut_save)
xline(45,LineWidth=4)
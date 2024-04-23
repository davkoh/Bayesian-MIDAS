%% The Local Trend Model with SV with Horseshoe Prior %%

clear all;

rng(1994,"twister");
% Choose Parameters of the DGP
    % Time-series
    T = 100;
    sd_tau = 0.3;
    t_0 = 0;
    sd_y = 1;
    h_0 = 0;
    sd_h = .1;
    % Group-correlated Data
    G = 5;
    K = 50;
    corr_x = 0.8;
    Sigma_X = kron(eye(G),ones(K/G)*corr_x);
    Sigma_X = Sigma_X- diag(diag(Sigma_X)) + diag(ones(K,1));
    X = mvnrnd(zeros(K,1),Sigma_X,T);
    beta = [1/G*ones(G,1);zeros(G,1);zeros(G,1);zeros(G,1);1/G*(ones(G,1));zeros(G,1);zeros(G,1);zeros(G,1);zeros(G,1);zeros(G,1)];
    beta_true = beta;


% Generate Data
tau = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);

outlier = 10;


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

nsim = 5000;
burnin = 5000;
valh = 0;
valg = 0;

%% prior
a0_h = 0; b0_h = 10;
a0_g = 0; b0_g = 10;
a0_tau = 0; b0_tau = 10;
Vomegah = .001;
Vomegag = .2;
    
% initialize the Markov chain
h0 = log(var(y))/5; g0 = log(var(y))/10; tau0 = mean(y);
omegah = sqrt(.2);
omegag = sqrt(.2);
h_tilde = zeros(T,1);
g_tilde = zeros(T,1);
h = h0 + omegah*h_tilde;
g = g0 + omegah*g_tilde; 
lamt = ones(T,1);
nut = 5;


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
beta_nu = .1;
nu_lb = 2;
nu_ub = 30; 
sig2_tau = .1; 
tau0 = 0;
a0 = 5; b0 = 100;
nu_sig0 = 3; S_sig0 = 1*(nu_sig0-1);
nu_omega0 = 3; S_omega0 = .25^2*(nu_omega0-1);

log_likelihood = @(x) -0.5*sum((y_ar - (2*x - 1).*x_ar).^2) + log(betapdf(x, prior_phi(1), prior_phi(2)));    
prior_phi = [10, 2]; % Parameters for Beta prior


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

   

%%
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


% define a few things
n_grid = 500; % number of grid points
omh_grid = linspace(-1,1,n_grid)';
omg_grid = linspace(-1,1,n_grid)';
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);

% initialize for storeage
store_theta = zeros(nsim,5); % [omegah omegag h0 g0 tau0]
store_tau = zeros(nsim,T); 
store_h = zeros(nsim,T);
store_g = zeros(nsim,T);
store_pomh = zeros(n_grid,1);
store_pomg = zeros(n_grid,1);
store_lpden = zeros(nsim,3); % log posterior densities of 
                             % omegah = valh; omegag = valg;
                             % omegah = omegag 
                            
rand('state', sum(100*clock) ); randn('state', sum(200*clock) );    
for isim = 1:nsim+burnin
    

    %% Sample beta
    data.Y = y-tau;
    data.iOh = 1./(exp(h).*lamt).*speye(T);
    [beta,hyper_params,data] = sample_beta(data,hyper_params);

        % sample tau    
    iOh = sparse(1:T,1:T,1./(exp(h).*lamt));
    HiOgH = H'*sparse(1:T,1:T,1./exp(g))*H;
    Ktau =  HiOgH + iOh;    
    tau_hat = Ktau\(tau0*HiOgH*ones(T,1) + iOh*(y-X*beta));
    tau = tau_hat + chol(Ktau,'lower')'\randn(T,1);
    
    % 2: Sample observation variance
        ystar = log((y-tau-X*beta).^2./lamt + 0.0001);
        [input_dyn_hs_sv] = sample_dyn_hs_sv(ystar,input_dyn_hs_sv);
         h = input_dyn_hs_sv.g;

        % sample g_tilde
    ystar = log((tau-[tau0;tau(1:end-1)]).^2 + .0001);
    [g_tilde g0 omegag omegag_hat Domegag] = ...
        SVRW_gam(ystar,g_tilde,g0,omegag,a0_g,b0_g,Vomegag); 
    g = g0 + omegag*g_tilde;
    
        % sample tau0
    Ktau0 = 1/b0_tau + 1/exp(g(1));
    tau0_hat = Ktau0\(a0_tau/b0_tau + tau(1)/exp(g(1)));
    tau0 = tau0_hat + sqrt(Ktau0)'\randn;

        %% Sample t-errors

    lamt = 1./gamrnd((nut+1)/2,2./(nut+(1./exp(h).*speye(T))*((y-tau-X*beta).^2)));

      nut = sample_nu_slice(nut,lamt,nu_lb,nu_ub,alpha,beta_nu);

            
    if (mod(isim, 1000) == 0)
        disp([num2str(isim) ' loops... ']);
    end     
    
    if isim > burnin
        isave = isim - burnin;
        store_tau(isave,:) = tau';
        store_h(isave,:) = h'; 
        store_g(isave,:) = g'; 
        store_theta(isave,:) = [omegah omegag h0 g0 tau0]; 
        save_nut_sw(isave) = nut;
            score_save(isave) = mean(normpdf(y,tau,lamt.*exp(h)));


        
    end    
end

theta_hat = mean(store_theta)';
tau_hat = mean(store_tau)';
h_hat = mean(exp(store_h/2))'; 
g_hat = mean(exp(store_g/2))'; 

clf;
plot(tau_hat)
hold on
plot(tau_true)

clf;
plot(exp(0.5*h_hat))
hold on
plot(exp(0.5*h_true))
plot(exp(0.5*g_true))
plot(exp(0.5*h_true) + exp(0.5*g_true))

trendrmse = sqrt(sum((tau_hat - tau_true).^2)/T)
svrmse = sqrt(sum((h_hat - h_true).^2)/T)



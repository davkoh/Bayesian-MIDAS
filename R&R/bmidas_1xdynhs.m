%% Function to do Inference with the HS on the Trend and SV %%


function out = bmidas_1xdynhs(input)

%% Unpack Data from Input Structure
% General Data
Y = input.Y;
X = input.X;
T = size(X,1);
K = size(X,2);
n_burn_in = input.burnin;
n_samples = input.samples;
% Model selection info
%trend_ind = input.trend;
%sv_ind = input.sv;
%t_ind =input.t;

% Data structure for the beta sampler

data.Y = input.Y;
data.X = input.X;
data.tX = data.X';
data.gl_param_expand_diag_inv = zeros(K,1);
data.K = K;
data.tau_shape_const = (K+1)/2;
data.stable_const = 1e-07;
data.G = size(unique(input.grp_idx),1);
data.grp_idx = input.grp_idx;
data.grp_size = histc(input.grp_idx, unique(input.grp_idx));
data.grp_size_cs = cumsum(histc(input.grp_idx, unique(input.grp_idx)));
data.a = input.a;
data.b = input.b;
data.iOh = eye(T);


%% Initialise Unknowns

% Prior for Trend SV
a0_g = 0; b0_g = 10;
a0_tau = 0; b0_tau = 10;
Vomegag = .2;
    
% initialize the Markov chain
g0 = log(var(Y))/10; tau0 = mean(Y);
omegag = sqrt(.2);
g_tilde = zeros(T,1);
g = g0 + omegag*g_tilde; 


    % For the beta sampler
    hyper_params.tau_sq = 1;
    hyper_params.nu = 5;
    hyper_params.gamma_sq = ones(size(unique(input.grp_idx),1),1);
    hyper_params.lambda_sq = ones(K,1);

    % Initalise series
tau = zeros(T,1);
phi = 0; % Start at zero perhaps
Hphi = speye(T) - sparse(2:T,1:(T-1),phi*ones(1,T-1),T,T);
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
mu = 0;
mu_tilde = [mu;(1-phi)*ones(T-1,1)*mu];
h = zeros(T,1);
h_tilde = h-mu;
xi = ones(T,1); % diagonals of the prior on h_tilde
iSig_xi = sparse(1:T,1:T,1./xi);
xi_0 = 0;
xi_mu = 0;
sig2_y = 1;
lamt = ones(T,1);
nut = 5;
alpha_nu = 2;
beta_nu = 1;
nu_lb = 2;
nu_ub = 30; 
% Function for log posterior of phi
log_likelihood = @(x) -0.5*sum((y_ar - (2*x - 1).*x_ar).^2) + log(betapdf(x, prior_phi(1), prior_phi(2)));    
prior_phi = [10, 2]; % Parameters for Beta prior
tau0 = 0;

V_a = 0.1;
V_b = 0.1;


% Prep input the sampling function
    input_dyn_hs_sv.g = h;
    input_dyn_hs_sv.phi = phi;
    input_dyn_hs_sv.g_tilde = h_tilde;
    input_dyn_hs_sv.Hphi = Hphi;
    input_dyn_hs_sv.iSig_xi = iSig_xi; 
    input_dyn_hs_sv.mu_tilde = mu_tilde;
    input_dyn_hs_sv.mu = mu;
    input_dyn_hs_sv.xi_mu = xi_mu;
    input_dyn_hs_sv.xi_0 = xi_0;
    input_dyn_hs_sv.y = Y;
    input_dyn_hs_sv.prior_phi = prior_phi; % Parameters for Beta prior
    input_dyn_hs_sv.H = H;


%% Storage Matrices
beta_store = zeros(K,n_samples);
gamma_store = zeros(size(unique(input.grp_idx),1),n_samples);
tausq_store = zeros(n_samples,1);
tau_store = zeros(T,n_samples);
store_h = zeros(T,n_samples);
store_g = zeros(T,n_samples);
nuy_store =zeros(n_samples,1);
lambda_store = zeros(T,n_samples);



for loops = 1:n_burn_in+n_samples
%% Sample beta
data.Y = Y-tau;
data.iOh = 1./(exp(h).*lamt).*speye(T);
[beta,hyper_params,data] = sample_beta(data,hyper_params);
%
        % sample tau    
    iOh = sparse(1:T,1:T,1./(exp(h).*lamt));
    HiOgH = H'*sparse(1:T,1:T,1./exp(g))*H;
    Ktau =  HiOgH + iOh;    
    tau_hat = Ktau\(tau0*HiOgH*ones(T,1) + iOh*(Y-X*beta));
    tau = tau_hat + chol(Ktau,'lower')'\randn(T,1);
    
    % 2: Sample observation variance
        ystar = log((Y-tau-X*beta).^2./lamt + 0.0001);
        [input_dyn_hs_sv] = sample_dyn_hs_sv(ystar,input_dyn_hs_sv);
         h = input_dyn_hs_sv.g;

        % sample g_tilde
    ystar = log((tau-[tau0;tau(1:end-1)]).^2 + .0001);
    [g_tilde g0 omegag omegag_hat Domegag] = ...
        SVRW_gam(ystar,g_tilde,g0,omegag,a0_g,b0_g,Vomegag); 
    %[g_tilde g0 omegag omegag_hat Domegag] = ...
    %    SVRW_gam_hier(ystar,g_tilde,g0,omegag,a0_g,b0_g,Vomegag,V_a,V_b); 
    g = g0 + omegag*g_tilde;
    
        % sample tau0
    Ktau0 = 1/b0_tau + 1/exp(g(1));
    tau0_hat = Ktau0\(a0_tau/b0_tau + tau(1)/exp(g(1)));
    tau0 = tau0_hat + sqrt(Ktau0)'\randn;

%% Sample t-errors
lamt = 1./gamrnd((nut+1)/2,2./(nut+(1./exp(input_dyn_hs_sv.g).*speye(T))*((Y-tau-X*beta).^2)));
%
nut = sample_nu_slice(nut,lamt,nu_lb,nu_ub,alpha_nu,beta_nu);
%[lamt,nut] = t_samp(Y,tau,X,beta,h,nut,nu_ub,H,lamt);
%
%% Save output

if loops>n_burn_in
    tau_store(:,loops-n_burn_in) = tau;
    tausq_store(loops-n_burn_in) = hyper_params.tau_sq;
    nuy_store(loops-n_burn_in) = nut;
    lambda_store(:,loops-n_burn_in) = lamt;
    beta_store(:,loops-n_burn_in) = beta;
    store_h(:,loops-n_burn_in) = input_dyn_hs_sv.g';
    store_g(:,loops-n_burn_in) = g';
    end
end

out.beta =beta_store;
out.tau = tau_store;
out.h = store_h;
out.g = store_g;
out.nu = nuy_store;
out.lambda = lambda_store;
out.beta_savs = savs(beta_store',X);



%% ============== Sampling Functions ============== %% 

% Put in the other functions here.


%% Save the SAVS sparsified Vectors
function [pbeta_savs] = savs(betaout,x)
for i = 1:size(betaout,1)
  for j = 1:size(x,2)
    xtx = x(:,j)'*x(:,j);
    mu = betaout(i,j)^(-2);
    if mu >= abs(betaout(i,j))*xtx
      pbeta_savs(i,j)=0;
    else
      pbeta_savs(i,j) = sign(betaout(i,j))*(abs(betaout(i,j))*xtx - mu)/xtx;
    end
  end
end

%% Sample t-distribution
function [lam,nu_y] =  t_samp(y,tau,X,beta,h,nu_y,nu_ub,H,lam) 

%% Description %%
% H = first difference matrix
% lam = inverse gamma distributed scales to force t-distribution
% nu_y = degrees of freedom of t-distribution

%% sample lam
temp1 = (H\(y-tau-X*beta)).^2./(exp(h).*lam)/2; 
lam =  1./gamrnd((nu_y+1)/2,2./(nu_y/2+temp1)); 

[nu_y,flag] = sample_nu_MH(lam,nu_y,nu_ub);

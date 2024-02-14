%% Function to do Inference with the HS on the Trend and SV %%


function out = bmidas_hspriors(input)

%% Unpack Data from Input Structure
% General Data
Y = input.Y;
X = input.X;
T = size(X,1);
K = size(X,2);
n_burn_in = input.burnin;
n_samples = input.samples;
% Model selection info
trend_ind = input.trend;
sv_ind = input.sv;
t_ind =input.t;

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

% Data structure for the trend sampler
prior_tau.V_t0 = 10;
prior_tau.H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);

% Data structure for the h sampler
prior_sv_h.Vh = 10;
% Data structure for the g sampler
prior_sv_g.Vh = 10;


%% Initialise Unknowns

    % For the beta sampler
    hyper_params.tau_sq = 1;
    hyper_params.nu = 5;
    hyper_params.gamma_sq = ones(size(unique(input.grp_idx),1),1);
    hyper_params.lambda_sq = ones(K,1);

    % For the Trend
    tau= zeros(T,1); % Trend
    Sigmatau = ones(T,1); % Covariance of trend
    tau0 = 0; % Initial condition
    lambdatau = ones(T,1); % Local shrinkage on trend
    nutau = 1; % Global shrinkage on trend

    g = ones(T,1);
    g0 = 0;
    Sigmag = ones(T,1);
    lambdag  = ones(T,1);
    nug = 1;
    a0_tau = 0; b0_tau = 10;

    % For Stochastic vol
    h = ones(T,1); % svol
    h0 = log(var(Y)); % initial condition svol
    Sigmah=ones(T,1); % covariance of vol
    lambdah=ones(T,1); % loval shrinkage on vol
    nuh=1; % global shrinkage on vol

    % For T-dist.
    nub = 50; % upper bound dof
    lamt = ones(T,1); % mixture for t-distribution
    nut = 5; % dof

    H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);


%% Storage Matrices
beta_store = zeros(K,n_samples);
gamma_store = zeros(size(unique(input.grp_idx),1),n_samples);
tausq_store = zeros(n_samples,1);
tau_store = zeros(T,n_samples);
store_h = zeros(T,n_samples);
nuy_store =zeros(n_samples,1);
lambda_store = zeros(T,n_samples);

%if trend_ind == 1
%store_alpha = zeros(n_samples,1); % Intercept if no trend is specified
%end

%% Adjustment to the Models Type

%if trend_ind == 0
%    Y = Y - mean(Y);
%end


%
for loops = 1:n_burn_in+n_samples
%% Sample beta
data.Y = Y;
data.iOh = 1./(exp(h).*lamt).*speye(T);
[beta,hyper_params,data] = sample_beta(data,hyper_params);

%% Sample Tau
iOh = sparse(1:T,1:T,1./(exp(h).*lamt));
HiOgH = H'*sparse(1:T,1:T,1./exp(g))*H;
Ktau =  HiOgH + iOh;    
tau_hat = Ktau\(tau0*HiOgH*ones(T,1) + iOh*(Y-tau));
tau = tau_hat + chol(Ktau,'lower')'\randn(T,1);

%% Sample g
ystar = log((tau-[tau0;tau(1:end-1)]).^2 + .0001);
[g,Sigmag,g0,lambdag,nug] = sample_sv_unitroot_hs(ystar,g,Sigmag,g0,lambdag,nug,prior_sv_g);

% Sample initial state:
    Ktau0 = 1/b0_tau + 1/exp(g(1));
    tau0_hat = Ktau0\(a0_tau/b0_tau + tau(1)/exp(g(1)));
    tau0 = tau0_hat + sqrt(Ktau0)'\randn;

%% Sample h
u = Y-tau-X*beta;
Ystar = log(u.^2.*lamt.^(-1) + .0001);
[h,Sigmah,h0,lambdah,nuh] = sample_sv_unitroot_hs(Ystar,h,Sigmah,h0,lambdah,nuh,prior_sv_h);

%% Sample t-errors
[lamt,nut] = sample_t_errors(u,lamt,nut,nub,exp(h));
lamt = ones(T,1);

%% Save output %% Needs adjustment to what is being saved

if loops>n_burn_in
    tau_store(:,loops-n_burn_in) = tau;
    tausq_store(loops-n_burn_in) = hyper_params.tau_sq;
    nuy_store(loops-n_burn_in) = nut;
    lambda_store(:,loops-n_burn_in) = lamt;
    beta_store(:,loops-n_burn_in) = beta;
    store_h(:,loops-n_burn_in) = h';
    end
end

out.beta =beta_store;
out.tau = tau_store;
out.h = store_h;
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
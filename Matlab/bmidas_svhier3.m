%% GGIG Regression Model Gibbs Sampler: Cleaned %%
% Rights reserved to David Kohns: david.kohns94@googlemail.com.
% X = TxM, covariate matrix on which to apply GIGG shrinkage
% C = TxK, covariate matrix on which to apply no shrinkage (intercept, other adjustment variables)
% Y = Tx1, response vector
% group_idx = 1xM vector which indicates which of the G groups the M
% covariates in X belong to
% grp_size = 1xG vector which indicates how large the individual groups are
% alpha_inits = Kx1 vector containing initial values for non-shrunk
% variables
% beta_inits = Mx1 vector containing initial values for regression
% coefficient vector
% lambda_sq_inits =
% gamma_sq_inits =
% p = Gx1 vector that contains the shape parameter for the GIG prior in the
% group parameters
% q= Gx1 vector that contains the shape parameter for the GIG prior in the
% inividual parameters
% tau_sq_init = 
% sigma_sq_init = 
% nu_init =
% n_burn_in = number of burnin samples
% n_samples = number of samples to save after burnin
% n_thin = 
% stable_const = parameter that controls numerical stability for the GIG
% posterior

function out = bmidas_svhier3(input)

warning('off','all')

addpath('/Users/dk/Documents/GitHub/MF-GIGG-Nowcasting-Project')

%% Unpack Data from Input Structure
grp_idx = input.grp_idx;
Y = input.Y;
X = input.X;
n_burn_in = input.burnin;
n_samples = input.samples;
a = input.a;
b = input.b;
standardise = input.standardise;
trend_ind = input.trend;
sv_ind = input.sv;
t_ind =input.t;

%% Precompute and store useful quantities
T = size(X,1);
G = size(unique(grp_idx),1); %% This might be difference between PC and Mac. Mac has spits out a 1 x G vector for unique(grp_idx)
K = size(X,2);
grp_size = histc(grp_idx, unique(grp_idx));
grp_size_cs = cumsum(grp_size);

if standardise == 1
    [X,mu_x,sig_x] = normalize(X);
    [Y,mu_y,sig_y] = normalize(Y);
end

tX = transpose(X);

% First difference matrix
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);




%% Initialise GIGG Stuff
beta = zeros(K,1);
lambda_sq = ones(K,1);
gamma_sq = ones(G,1);
tau_sq = 1;
nu = 5;
stable_const = 1e-07;

%% Initialise Trend-SV-t stuff

% Priors
a0_h = 0; b0_h = 10;
a0_g = 0; b0_g = 10;
a0_tau = 0; V_tau0 = 10;
Vomegah = 0.001;
Vomegag = 0.001;
nu_ub = 50;  % upper bound for nu
count_nu = 0;
V_a = .1;
V_b = .1;
Vh0 = .1;
Vg0 = .1;

% initialize the Markov chain
h0 = log(var(Y))/5; g0 = log(var(Y))/10; tau0 = mean(Y);
omegah = sqrt(.2);
omegag = sqrt(.2);
h_tilde = zeros(T,1);
g_tilde = zeros(T,1);
h = h0 + omegah*h_tilde;
g = g0 + omegah*g_tilde; 
nut = 5; %degrees of freedom
lamt = 1./gamrnd(nut/2,2/nut,T,1); % mixture weights for t
tau = zeros(T,1);
iOh = sparse(1:T,1:T,1./(exp(h).*lamt));

tau_shape_const = (K+1)/2;
gl_param_expand_diag_inv = zeros(K,1);

alpha_nu = 2;
beta_nu = 0.1;
nu_lb = 2;
nu_ub = 30; 
% Function for log posterior of phi
log_likelihood = @(x) -0.5*sum((y_ar - (2*x - 1).*x_ar).^2) + log(betapdf(x, prior_phi(1), prior_phi(2)));    
prior_phi = [10, 2]; % Parameters for Beta prior



%% Storage Matrices
beta_store = zeros(K,n_samples);
lambda_store = zeros(T,n_samples);
gamma_store = zeros(G,n_samples);
tau_store = zeros(T,n_samples);
tausq_store = zeros(n_samples,1);
sigma_store = zeros(n_samples,1);
nuy_store =zeros(n_samples,1);

store_theta = zeros(n_samples,8); % [omegah omegag h0 g0 tau0 V_omegah V_omegag V_tau0]
store_h = zeros(T,n_samples);
store_g = zeros(T,n_samples);
if trend_ind == 1
store_alpha = zeros(n_samples,1); % Intercept if no trend is specified
end

%% Adjustment to the Models Type

if trend_ind == 0
    Y = Y - mean(Y);
end


%
for loops = 1:n_burn_in+n_samples
%% Draw beta
[beta,tau_sq,gamma_sq,lambda_sq,gl_param_expand_diag_inv,nu] = beta_samp2(Y,tau,X,tX,iOh,tau_sq,gamma_sq,lambda_sq, ...
    gl_param_expand_diag_inv,K,tau_shape_const,stable_const,nu,G, ...
    grp_size_cs,grp_size,grp_idx,a,b, ...
    trend_ind,sv_ind,t_ind);
%}
%% Sample Tau

if trend_ind == 1

[tau] = trend_samp(Y,X,beta,iOh,g,H,tau0,T); %% Add a sv ind here as well
else
    tau = zeros(T,1);
end

%% Sample h_tilde

[h,iOh,h0,omegah,V_omegah,Vh0] = h_samp(Y,tau,X,beta,lamt,h_tilde,h0,omegah,a0_h,b0_h,Vomegah,Vh0,T,sv_ind,V_a,V_b);


if sv_ind == 1

%% Sample g_tilde

[g,g0,omegag, V_omegag,Vg0] = g_samp(tau,tau0,g_tilde,g0,omegag,a0_g,b0_g,Vomegag,Vg0,V_a,V_b);

%% Sample tau0
Ktau0 = 1/V_tau0 + 1/exp(g(1));
tau0_hat = Ktau0\(a0_tau/V_tau0 + tau(1)/exp(g(1)));
tau0 = tau0_hat + sqrt(Ktau0)'\randn;

% Sample Vtau0
V_tau0 = 1/gamrnd(V_a+1/2,1/(V_b+((tau0).^2)/2));

end

%% Sample parameters of the t-distribution 

if t_ind == 1

lamt = 1./gamrnd((nut+1)/2,2./(nut+(1./exp(input_dyn_hs_sv.g).*speye(T))*((Y-tau-X*beta).^2)));
%
nut = sample_nu_slice(nut,lamt,nu_lb,nu_ub,alpha_nu,beta_nu);
else
    lamt = ones(T,1);
    nut = 1;

end

%% Save output %% Needs adjustment to what is being saved

if loops>n_burn_in
    if trend_ind == 1
    tau_store(:,loops-n_burn_in) = tau;
    else
    store_alpha(loops-n_burn_in,:) = mean(input.Y) + sqrt(var(input.Y)/T)*randn;
    end
    tausq_store(loops-n_burn_in) = tau_sq;
    if t_ind == 1
    nuy_store(loops-n_burn_in) = nut;
    lambda_store(:,loops-n_burn_in) = lamt;
    end
    beta_store(:,loops-n_burn_in) = beta;
    %gamma_store(:,loops-n_burn_in) = gamma_sq;
    if sv_ind == 1
    store_h(:,loops-n_burn_in) = h'; 
    store_g(:,loops-n_burn_in) = g'; 
    store_theta(loops-n_burn_in,:) = [omegah omegag h0 g0 tau0 V_omegah V_omegag V_tau0]; 
    end
end

  %  if (mod(loops, 1000) == 0)
  %      disp([num2str(loops) ' loops of ' num2str(n_samples+n_burn_in)]);
  %  end 


end

out.beta =beta_store;
out.tau = tau_store;
out.h = store_h;
out.g = store_g;
out.theta = store_theta;
out.nu = nuy_store;
out.lambda = lambda_store;
out.beta_savs = savs(beta_store',X);
if trend_ind == 0
    out.alpha = store_alpha;
end
%}
%% ============== Sampling Functions ============== %% 

%% Sample Beta
function [beta,tau_sq,gamma_sq,lambda_sq,gl_param_expand_diag_inv,nu] = beta_samp2(Y,tau,X,tX,iOh,tau_sq,gamma_sq,lambda_sq, ...
    gl_param_expand_diag_inv,K,tau_shape_const,stable_const,nu,G, ...
    grp_size_cs,grp_size,grp_idx,a,b, ...
    trend_ind,sv_ind,t_ind)


if trend_ind == 1
yhat = Y; 
else
    yhat = Y;
end

if sv_ind == 1 || t_ind == 1
    iOh = iOh;
else
    iOh = eye(size(X,1));
end


    for gg  = 1:K
    gl_param_expand_diag_inv(gg) = 1.0 / (tau_sq * gamma_sq(grp_idx(gg)) * lambda_sq(gg));
    end
    beta_tmp = tX*iOh*X + sparse(diag(gl_param_expand_diag_inv)) + 1e-10;
    
    beta = tX*iOh*yhat + chol(beta_tmp,'lower')*randn(K,1);%beta = (1.0 / sigma_sq) * tX * (Y) + chol(beta_tmp,'lower')*randn(M,1);
   
    beta = beta_tmp\beta;


%% Draw tau^2
tau_rate_const = sum(beta.^2.*gl_param_expand_diag_inv);
tau_sq = 1.0 / gamrnd(tau_shape_const, 1.0 / (tau_sq * tau_rate_const / 2.0 + 1.0 / nu));

%% Draw gamma_g^2/lambda^2_gj
for j = 1:G

    % Sample gamma_g^2
    stable_psi = 0;
    if j == 1
    start_tmp = 1;
    end_tmp = grp_size_cs(j);
    else
    start_tmp = grp_size_cs(j-1)+1;
    end_tmp = grp_size_cs(j);
    end

    stable_psi = sum(beta(start_tmp:end_tmp).^2./lambda_sq(start_tmp:end_tmp));   
    stable_psi = stable_psi./tau_sq;
    stable_psi = max(stable_psi,stable_const);
    gamma_sq(j) = 1/gigrnd(grp_size(j)/2-a(j),stable_psi, 2, 1); %%%% Watch out for the p variable here.

    

    % Sample lambda^2_gj
    for i = 1:grp_size(j)
        lambda_sq(start_tmp+i-1) = 1.0 / gamrnd(b(j) + 0.5,...
            1.0 / (1 + (beta(start_tmp + i-1)^2) / (2.0 * tau_sq * gamma_sq(j))));
    end
end

%% Draw nu
nu = 1.0 / gamrnd(1, 1 / ((1 / tau_sq) ));

%% Sample Trend
function [tau] = trend_samp(y,X,beta,iOh,g,H,tau0,T)

%% Description %%
% iOh = covariance component due to SV-t in the observation equation
% tau0 = initial condition of the trend
% H = first difference matrix
%% Sample Tau
ystar = y-X*beta; 
HiOgH = H'*sparse(1:T,1:T,1./exp(g))*H;
Ktau =  HiOgH + iOh;    
tau_hat = Ktau\(tau0*HiOgH*ones(T,1) + iOh*ystar);
tau = tau_hat + chol(Ktau,'lower')'\randn(T,1);

%% Sample h
function [h,iOh,h0,omegah,Vomegah,Vh0] = h_samp(y,tau,X,beta,lam,h_tilde,h0,omegah,a0_h,b0_h,Vomegah,Vh0,T,sv_ind,V_a,V_b)

%% Description %%
% iOh = covariance component due to SV-t in the observation equation
% tau0 = initial condition of the trend
% H = first difference matrix
%% Sample h

if sv_ind ==1
ystar =  log(((y-tau-X*beta)).^2 + .0001); 
[h_tilde h0 omegah omegah_hat Domegah Vomegah Vh0] = ...
    SVRW_gam_hier2(ystar,h_tilde,h0,omegah,a0_h,b0_h,Vomegah,Vh0,V_a,V_b); 
h = h0 + omegah*h_tilde;   

else
    h = ones(T,1);
end

iOh = sparse(1:T,1:T,1./(exp(h).*lam));

%% Sample g
function [g,g0,omegag,Vomegag,Vg0] = g_samp(tau,tau0,g_tilde,g0,omegag,a0_g,b0_g,Vomegag,Vg0,V_a,V_b)

ystar = log((tau-[tau0;tau(1:end-1)]).^2 + .0001);
[g_tilde g0 omegag omegag_hat Domegag Vomegag Vg0] = ...
    SVRW_gam_hier2(ystar,g_tilde,g0,omegag,a0_g,b0_g,Vomegag,Vg0,V_a,V_b); 
g = g0 + omegag*g_tilde;


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



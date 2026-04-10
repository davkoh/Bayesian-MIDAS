%% Horseshoe Regression Model Gibbs Sampler

function out = bmidas_horseshoe(input)

warning('off','all')


%% Unpack Data from Input Structure

% Model Structure
    % Which MIDAS
midas_prior = input.prior.midas;
    % SV-Obs type
sv_obs_type = input.prior.sv_obs;
    % Which trend type
trend_type =  input.prior.trend_sv;
    % Which tail type 
tail_type = input.prior.tail_type;

% Prior Structure
gigg_type =  input.prior.gigg_hyper;
ag = input.prior.a_g;
bg = input.prior.b_g;
V_omegag = input.prior.V_omegag;
V_g0 = input.prior.V_g0;
V_tau0 = input.prior.V_tau0;
xi_g = input.prior.xi_g; % where is that supposed to be used? 
V_omegah =  input.prior.V_omegah;
V_h0 = input.prior.V_h0;
xi_h = input.prior.xi_h ; % where is that supposed to be used?

% Data for model
grp_idx = input.grp_idx;
Y = input.Y;
X = input.X;
standardise = input.standardise;
if strcmp(trend_type,"none") == 1
    Y = Y - mean(Y);
end


% MCMC Sampler
n_burn_in = input.burnin;
n_samples = input.samples;


% Further priors assumed to be unchanged
    % For t-distribution t_nu(exp(h_t))
nu_ub = 30;  % upper bound for nu
nu_lb = 2; % lower bound for nu (2 needed for 2 finite fractional moments)
rate1_t = 2; % recommendation of the paper
rate2_t = 0.1; % recommendation of the paper

    % For GIGG(a_g,b_g) hyperparameters if hierarchical GIGG is selected
    % (recommendation of the paper)
rate1 = 1;
rate2 = 3;



%% Precompute and store useful quantities
T = size(X,1);
G = size(unique(grp_idx),1);
K = size(X,2);
grp_size = histc(grp_idx, unique(grp_idx));
grp_size_cs = cumsum(grp_size);

if standardise == 1
    [X,mu_x,sig_x] = normalize(X);
    [Y,mu_y,sig_y] = normalize(Y);
end

tX = transpose(X);



%% Initialise chains for MCMC

% Parameters related to the MIDAS component 
theta = zeros(K,1); % MIDAS coefficients
lambdabeta=ones(K,1);
taubeta = 1;
nubeta = ones(K,1);
etabeta = 1;
iVbeta = diag(ones(K,1));

% Parameters related to the trend
tau = zeros(T,1);
h0 = log(var(Y))/5; g0 = log(var(Y))/10; tau0 = mean(Y);
omegah = sqrt(.2);
omegag = sqrt(.2);
h_tilde = zeros(T,1);
g_tilde = zeros(T,1);

if strcmp(sv_obs_type,"none") 
    h = zeros(T,1);
    g = zeros(T,1);
else
    h = h0 + omegah*h_tilde;
    g = g0 + omegah*g_tilde;

end

% Parameters related to the t-distribution
nu_y = 6; %degrees of freedom of t-distribution

if  strcmp(tail_type,"terr") 
lam = 1./gamrnd(nu_y/2,2/nu_y,T,1); % mixture weights for t
else
    lam = ones(T,1);
end

% Define joint covariance of observation error
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T); % difference matrix

if ~strcmp(sv_obs_type,"none")
iOh = sparse(1:T,1:T,1./(exp(h).*lam));
else
    iOh = sparse(1:T,1:T,1);
end

%% Storage Matrices

theta_store = zeros(K,n_samples); % MIDAS coefficients
varphi_store = zeros(T,n_samples); % Local variances
tau_store = zeros(T,n_samples); % Trend 
varthetasq_store = zeros(n_samples,1); % Global variances
sigma_store = zeros(n_samples,1); % observation variance, when not SV is selected
nuy_store =zeros(n_samples,1); % degrees of freedom of the normal

store_ag = zeros(G,n_samples); % hyper-parameter for gamma_k
store_bg = zeros(G,n_samples); % hyper-parameter for vartheta_{k,j}
store_state_params = zeros(n_samples,5); % [omegah omegag h0 g0 tau0]
store_h = zeros(T,n_samples); % stochastic volatility for the observation equation
store_g = zeros(T,n_samples); % stochastic volatility for the trend equation

if strcmp(trend_type,"none") ==1
store_alpha = zeros(n_samples,1); % Intercept if no trend is specified
end


for loops = 1:n_burn_in+n_samples
   
   yhat = Y; 

   if sv_ind == 1 || t_ind == 1
    iOh = iOh;
else
    iOh = sparse(1:T,1:T,1./(sigma_sq));
   end
    
    Kbeta = sparse(iVbeta) + X'*iOh*X  ; 
    theta_hat = Kbeta\(X'*iOh*yhat);
    theta = theta_hat + chol(Kbeta,'lower')'\randn(size(X,2),1);

    % sample lambdabeta
    lambdabeta = 1./gamrnd(1, 1./( 1./nubeta + 0.5*theta.^2/taubeta ));
    % sample taubeta
    taubeta = 1/gamrnd( 0.5*(size(X,2)), 1/( 1/etabeta + 0.5*sum(sum(beta.^2./lambdabeta))  ) );
    % sample nubeta
    nubeta = 1./gamrnd(1, 1./(1 + 1./lambdabeta));
    % samplel etabeta
    etabeta = 1/gamrnd(1, 1/( 1 + 1/taubeta ));
    iVbeta = diag(taubeta*lambdabeta)\speye(size(X,2));
    


    % Draw vartheta^2
     tau_rate_const = sum(theta.^2.*gl_param_expand_diag_inv);
     vartheta_sq = 1.0 / gamrnd(tau_shape_const, 1.0 / (vartheta_sq * tau_rate_const / 2.0 + 1.0 / nu));


%% Sample trend

if strcmp(trend_type,"fixed_SV") ==1 || strcmp(trend_type,"PC") == 1

y_star = Y-X*theta;
HiOgH = H'*sparse(1:T,1:T,1./exp(g))*H;
Ktau =  HiOgH + iOh;    
tau_hat = Ktau\(tau0*HiOgH*ones(T,1) + iOh*y_star);
tau = tau_hat + chol(Ktau,'lower')'\randn(T,1);

else

    tau = zeros(T,1);
end

%% Sample h_tilde:


ystar = log((Y-tau-X*theta).^2./lam + .0001);

% Fixed SV observation equation
if strcmp(sv_obs_type,"fixed_SV") == 1

[h_tilde h0 omegah omegah_hat Domegah] = ...
    SVRW_gam_omori(ystar,h_tilde,h0,omegah,0,V_h0,V_omegah); 
h = h0 + omegah*h_tilde;  

end

% PC SV observation equation
if strcmp(sv_obs_type,"PC") == 1

[h_tilde,h0,omegah,V_omegah,V_h0] = ...
    SVRW_gam_omori_pc(ystar,h_tilde,h0,omegah,0,V_omegah,V_h0);
h = h0 + omegah*h_tilde;  

end


if strcmp(sv_obs_type,"none")
    sigma_sq = 1/gamrnd((T+1)/2,1/((Y  - X * theta-tau)'*(Y - X * theta-tau)/2 ));
    h = zeros(T,1); %% Changed to zero if none
end

% Update var-covar of the observation equation
if strcmp(sv_obs_type,"none")
    iOh = sparse(1:T,1:T,1./(sigma_sq));
else
   iOh = sparse(1:T,1:T,1./(exp(h).*lam));
end




%% Sample g_tilde: TODO: add PC prior function
ystar = log((tau-[tau0;tau(1:end-1)]).^2 + .0001);

if strcmp(trend_type,"fixed_SV") ==1
    ystar = log((tau-[tau0;tau(1:end-1)]).^2 + .0001);
    [g_tilde g0 omegag omegag_hat Domegag] = ...
    SVRW_gam_omori(ystar,g_tilde,g0,omegag,0,V_g0,V_omegag); 

    g = g0 + omegag*g_tilde;
end

if strcmp(trend_type,"PC") ==1
    
    [g_tilde,g0,omegag,V_omegag,V_g0] = ...
    SVRW_gam_omori_pc(ystar,g_tilde,g0,omegag,0,V_omegag,V_g0);

    g = g0 + omegag*g_tilde;  

end

if strcmp(trend_type,"PC") ==1 || strcmp(trend_type,"fixed_SV") ==1
% Sample tau0
Ktau0 = 1/V_tau0 + 1/exp(g(1));
tau0_hat = Ktau0\(0/V_tau0 + tau(1)/exp(g(1)));
tau0 = tau0_hat + sqrt(Ktau0)'\randn;
end


if strcmp(trend_type,"none")
    g = zeros(T,1);
end


%% sample t-distribution parameters

if strcmp(tail_type,"terr")
e = Y - X*theta -tau;
lam = 1./gamrnd((nu_y+1)/2,2./(nu_y+e.^2./exp(h)));
%
nu_y = sample_nu_slice(nu_y,lam,nu_lb,nu_ub,rate1_t,rate2_t);
end


%% Save output

if loops>n_burn_in
    if strcmp(trend_type,"fixed_SV") == 1 || strcmp(trend_type,"PC")  == 1
    tau_store(:,loops-n_burn_in) = tau;
    store_g(:,loops-n_burn_in) = g'; 
    else
    store_alpha(loops-n_burn_in,:) = mean(input.Y) + sqrt(var(input.Y)/T)*randn; % non-informative prior for the intercept
    end
    varthetasq_store(loops-n_burn_in) = vartheta_sq;
    if strcmp(tail_type,"terr")
    nuy_store(loops-n_burn_in) = nu_y;
    %varphi_store(:,loops-n_burn_in) = varphi_sq;
    end
    theta_store(:,loops-n_burn_in) = theta;
    if ~strcmp(sv_obs_type,"none")
    store_h(:,loops-n_burn_in) = h'; 
    store_state_params(loops-n_burn_in,:) = [omegah omegag h0 g0 tau0]; 
    else
        sigma_store(loops-n_burn_in) = sigma_sq;
    end
    store_ag(:,loops-n_burn_in) = ag;
    store_bg(:,loops-n_burn_in) = bg;
end

    if (mod(loops, 1000) == 0)
        disp([num2str(loops) ' loops... ']);
    end 


end

out.theta =theta_store;
out.tau = tau_store;
out.h = store_h;
out.g = store_g;
out.state_params = store_state_params;
out.nu = nuy_store;
out.varphi = varphi_store;
out.sigma2 = sigma_store;
if strcmp(trend_type,"none") == 1 
    out.alpha = store_alpha;
end
out.a = store_ag;
out.b = store_bg;
%out.tau_varinv = squeeze(mean(store_ktauinv,3));
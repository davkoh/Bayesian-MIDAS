%% Horseshoe Regression Model Gibbs Sampler

function out = bmidas_mal(input)

warning('off','all')


%% Unpack Data from Input Structure



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
varphi_sq = ones(K,1);
gamma_sq = ones(G,1);
vartheta_sq = 1;
sigma_sq = var(Y);
nu = 1; % mixture variable for the Cauchy distribution of the global variance parameter in the GIGG (vartheta_sq)
stable_const = 1e-07;
stable_hyp_lb = 0.0001; % lower bound for ag posterior (Boss et al. (2024))
stable_hyp_ub = 1; % upper bound for bg posterior (Boss et al. (2024))
tau_shape_const = (K+1)/2; 
gl_param_expand_diag_inv = zeros(K,1);

% Parameters related to the trend
tau = zeros(T,1);
h0 = log(var(Y))/5; g0 = log(var(Y))/10; tau0 = mean(Y);
omegah = sqrt(.2);
omegag = sqrt(.2);
h_tilde = zeros(T,1);
g_tilde = zeros(T,1);

if sv_ind == 1 % TODO: check whether sv_ind is still needed, with the sv_type variable
h = h0 + omegah*h_tilde;
g = g0 + omegah*g_tilde;
else
    h = ones(T,1);
    g = ones(T,1);

end

% Parameters related to the t-distribution
nu_y = 6; %degrees of freedom of t-distribution

if  t_ind ==1
lam = 1./gamrnd(nu/2,2/nu,T,1); % mixture weights for t
else
    lam = ones(T,1);
end

% Define joint covariance of observation error
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T); % difference matrix

if sv_ind == 1 || t_ind ==1
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
store_ktauinv = zeros(T,T,n_samples); % covar of the trend TODO: is this still needed? 

store_ag = zeros(G,n_samples); % hyper-parameter for gamma_k
store_bg = zeros(G,n_samples); % hyper-parameter for vartheta_{k,j}
store_state_params = zeros(n_samples,5); % [omegah omegag h0 g0 tau0]
store_h = zeros(T,n_samples); % stochastic volatility for the observation equation
store_g = zeros(T,n_samples); % stochastic volatility for the trend equation

if strcmp(trend_type,"none") ==1
store_alpha = zeros(n_samples,1); % Intercept if no trend is specified
end


for loops = 1:n_burn_in+n_samples
    %% Draw theta (MIDAS coefficients)
   yhat = Y; 

   if sv_ind == 1 || t_ind == 1
    iOh = iOh;
else
    iOh = sparse(1:T,1:T,1./(sigma_sq));
   end

    for gg  = 1:K
    gl_param_expand_diag_inv(gg) = 1.0 / (vartheta_sq * gamma_sq(grp_idx(gg)) * varphi_sq(gg));
    end
    theta_tmp = tX*iOh*X + sparse(diag(gl_param_expand_diag_inv)) + 1e-10;
    theta = tX*iOh*yhat + chol(theta_tmp,'lower')*randn(K,1);%theta = (1.0 / sigma_sq) * tX * (Y) + chol(theta_tmp,'lower')*randn(M,1);
    theta = theta_tmp\theta;


    % Draw vartheta^2
        tau_rate_const = sum(theta.^2.*gl_param_expand_diag_inv);
        vartheta_sq = 1.0 / gamrnd(tau_shape_const, 1.0 / (vartheta_sq * tau_rate_const / 2.0 + 1.0 / nu));


    % Draw gamma_k^2/varphi^2_kj:
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
    
    % Option for inference on hyper-parameter, a_g
   if  strcmp(gigg_type,"hier_ag")==1 || strcmp(gigg_type,"hier_ag_bg")==1  
    ag(j) = sample_ag_slice(ag(j),gamma_sq(j),stable_hyp_lb,stable_hyp_ub,rate1,rate2); % Samples the hierarchical a_g component 
    end

    stable_psi = sum(theta(start_tmp:end_tmp).^2./varphi_sq(start_tmp:end_tmp));   
    stable_psi = stable_psi./vartheta_sq;
    stable_psi = max(stable_psi,stable_const);
    gamma_sq(j) = 1/gigrnd(grp_size(j)/2-ag(j),stable_psi, 2, 1); 


    % Sample varphi^2_kj
    for i = 1:grp_size(j)
        varphi_sq(start_tmp+i-1) = 1.0 / gamrnd(bg(j) + 0.5,...
            1.0 / (1 + (theta(start_tmp + i-1)^2) / (2.0 * vartheta_sq * gamma_sq(j))));
    end

    % Option for inference on hyper-parameter, b_g
    if strcmp(gigg_type,"hier_bg")==1 || strcmp(gigg_type,"hier_ag_bg")==1  
    bg(j) = sample_bg_slice(bg(j),varphi_sq(start_tmp:end_tmp),stable_hyp_lb,stable_hyp_ub,rate1,rate2);
    end
end

% Draw nu (mixture variable for Cauchy distribution)
nu = 1.0 / gamrnd(1, 1 / ((1 / vartheta_sq) ));

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

%% Sample h_tilde: TODO: add PC prior function


ystar = log((Y-tau-X*theta).^2./lam + .0001);

% Fixed SV-trend
if strcmp(trend_type,"fixed_SV") == 1

[h_tilde h0 omegah omegah_hat Domegah] = ...
    SVRW_gam_omori(ystar,h_tilde,h0,omegah,0,V_h0,V_omegah); 
h = h0 + omegah*h_tilde;  

end

% PC SV-trend. 
if strcmp(trend_type,"PC") == 1

[h_tilde,h0,omegah,V_omegah,V_h0] = ...
    SVRW_gam_omori_pc(ystar,h_tilde,h0,omegah,0,V_omegah,V_h0);
h = h0 + omegah*h_tilde;  

end


if sv_ind ==0
    sigma_sq = 1/gamrnd((T+1)/2,1/((Y  - X * theta-tau)'*(Y - X * theta-tau)/2 ));
    h = ones(T,1);
end

% Update var-covar of the observation equation
if sv_ind || t_ind == 1

iOh = sparse(1:T,1:T,1./(exp(h).*lam));
else
   iOh = sparse(1:T,1:T,1./(sigma_sq));
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
    g = ones(T,1);
end


%% sample t-distribution parameters

if t_ind ==1
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
    if t_ind == 1
    nuy_store(loops-n_burn_in) = nu_y;
    varphi_store(:,loops-n_burn_in) = varphi_sq;
    end
    theta_store(:,loops-n_burn_in) = theta;
    if sv_ind == 1
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
if strcmp(trend_type,"fixed") == 1 
    out.alpha = store_alpha;
end
out.a = store_ag;
out.b = store_bg;
%out.tau_varinv = squeeze(mean(store_ktauinv,3));
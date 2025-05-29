%% GGIG Regression Model Gibbs Sampler


function out = gigg_bsts_sv_t_aut(input)

warning('off','all')


%% Unpack Data from Input Structure
grp_idx = input.grp_idx;
Y = input.Y;
X = input.X;
n_burn_in = input.burnin;
n_samples = input.samples;
p = input.a;
q = input.b;
standardise = input.standardise;
sv_ind = input.sv_ind;
t_ind = input.t_ind;
trend_ind = input.trend_ind;




%% Precompute and store useful quantities
T = size(X,1);
G = size(unique(grp_idx),1);
M = size(X,2);
K = size(X,2);
grp_size = histc(grp_idx, unique(grp_idx));
grp_size_cs = cumsum(grp_size);

if standardise == 1
    [X,mu_x,sig_x] = normalize(X);
    [Y,mu_y,sig_y] = normalize(Y);
end

tX = transpose(X);



%% Initialise GIGG Stuff
beta = zeros(M,1);
lambda_sq = ones(M,1);
gamma_sq = ones(G,1);
tau_sq = 1;
sigma_sq = var(Y);
nu = 1;
stable_const = 1e-07;

%% Initialise Trend-SV-t stuff

% Priors
a0_h = 0; b0_h = 10;
a0_g = 0; b0_g = 10;
a0_tau = 0; b0_tau = 10;
Vomegah = .001;
Vomegag = .001;
nu_ub = 30;  % upper bound for nu

% initialize the Markov chain
h0 = log(var(Y))/5; g0 = log(var(Y))/10; tau0 = mean(Y);
omegah = sqrt(.2);
omegag = sqrt(.2);
h_tilde = zeros(T,1);
g_tilde = zeros(T,1);

if sv_ind == 1 
h = h0 + omegah*h_tilde;
g = g0 + omegah*g_tilde;
else
    h = ones(T,1);
    g = ones(T,1);

end

nu_y = 6; %degrees of freedom

if t_ind == 1 
lam = 1./gamrnd(nu/2,2/nu,T,1); % mixture weights for t
else
    lam = ones(T,1);
end

% define a few things

tau = zeros(T,1);

if sv_ind == 1 || t_ind ==1
iOh = sparse(1:T,1:T,1./(exp(h).*lam));
else
    iOh = sparse(1:T,1:T,1);
end

%% Storage Matrices
beta_store = zeros(K,n_samples);
lambda_store = zeros(T,n_samples);
tau_store = zeros(T,n_samples);
tausq_store = zeros(n_samples,1);
sigma_store = zeros(n_samples,1);
nuy_store =zeros(n_samples,1);

store_theta = zeros(n_samples,5); % [omegah omegag h0 g0 tau0]
store_h = zeros(T,n_samples);
store_g = zeros(T,n_samples);
if trend_ind == 0
store_alpha = zeros(n_samples,1); % Intercept if no trend is specified
end




%% Calculate constants for updating sigma and tau
tau_shape_const = (M+1)/2;

%% Prevent repatative initialisations by initialising here

gl_param_expand_diag_inv = zeros(M,1);
 

if trend_ind == 0
    Y = Y - mean(Y);
end

for loops = 1:n_burn_in+n_samples
    %% Draw beta
   yhat = Y; 

   if sv_ind == 1 || t_ind == 1
    iOh = iOh;
else
    iOh = eye(size(X,1));
   end

   %% Sample beta

    for gg  = 1:M
    gl_param_expand_diag_inv(gg) = 1.0 / (tau_sq * gamma_sq(grp_idx(gg)) * lambda_sq(gg));
    end
    beta_tmp = tX*iOh*X + sparse(diag(gl_param_expand_diag_inv)) + 1e-10;
    beta = tX*iOh*yhat + chol(beta_tmp,'lower')*randn(M,1);
    beta = beta_tmp\beta;


% Draw tau^2
tau_rate_const = sum(beta.^2.*gl_param_expand_diag_inv);
tau_sq = 1.0 / gamrnd(tau_shape_const, 1.0 / (tau_sq * tau_rate_const / 2.0 + 1.0 / nu));


% Draw gamma_g^2/lambda^2_gj: Just use the gig function in the matfiles.
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
    gamma_sq(j) = 1/gigrnd(grp_size(j)/2-p(j),stable_psi, 2, 1); %%%% Watch out for the p variable here.


    % Sample lambda^2_gj
    for i = 1:grp_size(j)
        lambda_sq(start_tmp+i-1) = 1.0 / gamrnd(q(j) + 0.5,...
            1.0 / (1 + (beta(start_tmp + i-1)^2) / (2.0 * tau_sq * gamma_sq(j))));
    end
end

% Draw nu
nu = 1.0 / gamrnd(1, 1 / ((1 / tau_sq) ));

%% Sample trend

if trend_ind == 1


y_star = Y-X*beta;

HiOgH = H'*sparse(1:T,1:T,1./exp(g))*H;
Ktau =  HiOgH + iOh;    
tau_hat = Ktau\(tau0*HiOgH*ones(T,1) + iOh*y_star);
tau = tau_hat + chol(Ktau,'lower')'\randn(T,1);

else

    tau = zeros(T,1);
end

%% Sample h_tilde


ystar = log((Y-tau-X*beta).^2./lam + .0001);
[h_tilde h0 omegah omegah_hat Domegah] = ...
    SVRW_gam_omori(ystar,h_tilde,h0,omegah,a0_h,b0_h,Vomegah); 

if sv_ind ==0
    sigma_sq = 1/gamrnd((T+1)/2,1/((Y  - X * beta-tau)'*(Y - X * beta-tau)/2 ));
end

if sv_ind == 1

h = h0 + omegah*h_tilde;   

else
    h = ones(T,1);
end

if sv_ind || t_ind == 1

iOh = sparse(1:T,1:T,1./(exp(h).*lam));
else
    sparse(1:T,1:T,1./(sigma_sq));
end



%% Sample g_tilde

if trend_ind == 1
ystar = log((tau-[tau0;tau(1:end-1)]).^2 + .0001);
[g_tilde g0 omegag omegag_hat Domegag] = ...
    SVRW_gam_omori(ystar,g_tilde,g0,omegag,a0_g,b0_g,Vomegag); 

g = g0 + omegag*g_tilde;

% Sample tau0
Ktau0 = 1/b0_tau + 1/exp(g(1));
tau0_hat = Ktau0\(a0_tau/b0_tau + tau(1)/exp(g(1)));
tau0 = tau0_hat + sqrt(Ktau0)'\randn;

else
    g = ones(T,1);
end


%% sample lam
if t_ind ==1
e = Y - X*beta -tau;
lam = 1./gamrnd((nu_y+1)/2,2./(nu_y+e.^2./exp(h)));
%
nu_y = sample_nu_slice(nu_y,lam,2,nu_ub,2,1);
end


%% Save output

if loops>n_burn_in
    if trend_ind == 1
    tau_store(:,loops-n_burn_in) = tau;
    else
    store_alpha(loops-n_burn_in,:) = mean(input.Y) + sqrt(var(input.Y)/T)*randn;
    end
    tausq_store(loops-n_burn_in) = tau_sq;
    if t_ind == 1
    nuy_store(loops-n_burn_in) = nu_y;
    lambda_store(:,loops-n_burn_in) = lam;
    end
    beta_store(:,loops-n_burn_in) = beta;
    if sv_ind == 1
    store_h(:,loops-n_burn_in) = h'; 
    store_g(:,loops-n_burn_in) = g'; 
    store_theta(loops-n_burn_in,:) = [omegah omegag h0 g0 tau0]; 
    else
        sigma_store(loops-n_burn_in) = sigma_sq;
    end
end

    if (mod(loops, 1000) == 0)
        disp([num2str(loops) ' loops... ']);
    end 


end

out.beta =beta_store;
out.tau = tau_store;
out.h = store_h;
out.g = store_g;
out.theta = store_theta;
out.nu = nuy_store;
out.lambda = lambda_store;
out.sigma2 = sigma_store;
if trend_ind == 0
    out.alpha = store_alpha;
end
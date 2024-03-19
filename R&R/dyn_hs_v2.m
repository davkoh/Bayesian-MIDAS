%% The Local Trend Model with SV with Horseshoe Prior %%

clear all;

% Choose Parameters of the DGP
T = 200;
sd_tau = 0.01;
t_0 = 0;
sd_y = 0.1;

% Generate Data
tau = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);

t = [0.1:0.1:20]
a = sin(t);

for t= 1:T
    if t == 1
        tau(1) = t_0 + sd_tau*randn;
    else
        tau(t) = a(t-1) + sd_tau*randn;
    end

    if t > 100
        tau(t) = a(t-1) + sd_tau*randn;
        y(t) = tau(t) + 10 + sd_y*randn;
    else
        y(t) = tau(t) + sd_y*randn;
    end
end


tau_true = tau;
h_true = h;

clf;
plot(y)
plot(tau_true)

% Priors
nu_sig0 = 3; S_sig0 = 1*(nu_sig0-1);


% Initalise series
tau0 = 0;
lam2t = zeros(T,1);
nu2 = 0;
phi = 0.9;
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

% Sampler Specs
rng(1,'twister');
iter = 20000;
burn = iter/2;
mcmc = iter-burn;

% Function for log posterior of phi
log_likelihood = @(x) -0.5*sum((y_ar - (2*x - 1).*x_ar).^2) + log(betapdf(x, prior_phi(1), prior_phi(2)));    
prior_phi = [10, 2]; % Parameters for Beta prior


%% Gibbs Sampler

% Storage
tau_save = NaN(T,mcmc);
tau0_save = NaN(mcmc,1);
sig2_save = NaN(mcmc,1);


for i = 1:iter
    %% 1: Sample Dynamic Shrinkage Components
        % (a) Sample h_tilde

        % Sample h_t+1 = mu + phi*(h_t-mu) + eta_t
        ystar = log([tau(1)-tau0;tau(2:end)-tau(1:end-1)].^2 + 0.0001);

        % normal mixture
        pj = [.22715 .20674 .18842 .13057 .12047 .05591 .04775 .01575	.00609, .00115];
        mj = [-0.85173 .02266 -1.97278 .73504 -3.46788 -5.55246 1.34744 -8.68384 1.92677 -14.65];  % warning: means already adjusted
        sigj2 = [.62699 .40611 .98583 .26768 1.57469 2.54498 .17788 4.16591 .11265 7.33342];
        sigj = sqrt(sigj2);
        % sample S from a 10-point distrete distribution
        temprand = rand(T,1);
        q = repmat(pj,T,1).*normpdf(repmat(ystar,1,10),...
        repmat(h_tilde,1,10)+repmat(mj,T,1),repmat(sigj,T,1));
        q = q./repmat(sum(q,2),1,10);
        S = 10 - sum(repmat(temprand,1,10)<cumsum(q,2),2)+1;
        m_s = mj(S)'; iSig_v = sparse(1:T,1:T,1./sigj2(S));

        Kh= (iSig_v + Hphi'*iSig_xi*Hphi);
        Ch = chol(Kh,'lower');
        htilde_hat = Kh\(iSig_v*(ystar- m_s -mu_tilde));
        h_tilde = htilde_hat + Ch'\randn(T,1);
        h = h_tilde + mu;

        
        % (b) Sample xi_t (PG variables)
        eta = Hphi*h_tilde;
        xi = pgdraw(eta);
        iSig_xi = sparse(1:T,1:T,1./xi);

        % (c) Sample mu
        Kmu = xi_mu + xi_0 + (1-phi^2)*sum(xi(1:T-1));
        Cm = chol(Kmu,'lower');
        mu_hat = Kmu\((xi_0*h(1) + (1-phi)*sum(xi(1:T-1).*[h(2:T)-phi*h(1:T-1)])));
        mu = mu_hat + Cm'\randn(1,1);
        mu_tilde = [mu;(1-phi)*ones(T-1,1)*mu];

        % (d) Sample xi_mu
        xi_mu = pgdraw(mu);

        % (e) Sample xi_0 
        xi_0 = pgdraw(h_tilde(1));

        % (f) Sample phi (needs slice sampling)
            % define y_ar
            y_ar = h_tilde(2:end)./(xi(2:end).^(1/2));
            x_ar = h_tilde(1:end-1)./(xi(1:end-1).^(1/2));

            % define x_ar
        phi = slice_sampling(phi, log_likelihood, 0, 1); % Adjust upper and lower bounds as needed
        Hphi = speye(T) - sparse(2:T,1:(T-1),phi*ones(1,T-1),T,T);


    %% 2: Sample State variables
        % (a) tau
        iSig_tau = sparse(1:T,1:T,1./exp(h));
        Ktau= (eye(T)*1/sig2_y + H'*iSig_tau*H);
        Ctau = chol(Ktau,'lower');
        tau_hat = Ktau\(1/sig2_y*y);
        tau = tau_hat + Ctau'\randn(T,1);

        % (b) tau_0 (not needed)

    

    %% 3: Sample observation variance
        sig2_y = 1/gamrnd(nu_sig0 + T/2,1/(S_sig0 + (y-tau)'*(y-tau)/2));
       

if i> burn 
    isave = i - burn;
    tau_save(:,isave) = tau;
    tau0_save(isave) = phi;
    
end

end

ptau = mean(tau_save,2);
clf;
plot(ptau)
hold on
plot(y)
plot(y)

sqrt(sum((ptau-y).^2)/T)

ph = mean(h_save,2);
clf;
plot(exp(0.5*ph))
hold on
plot(exp(0.5*h_true))

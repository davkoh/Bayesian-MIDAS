function [output_dyn_hs_trend] = sample_dyn_hs_trend(ystar,input_dyn_hs_trend)

% Define useful things
T = length(ystar);

% Unpack input
    tau = input_dyn_hs_trend.tau;
    h = input_dyn_hs_trend.h;
    phi = input_dyn_hs_trend.phi;
    h_tilde = input_dyn_hs_trend.h_tilde;
    Hphi = input_dyn_hs_trend.Hphi;
    iSig_xi = input_dyn_hs_trend.iSig_xi; 
    mu_tilde  = input_dyn_hs_trend.mu_tilde;
    mu = input_dyn_hs_trend.mu;
    xi_mu = input_dyn_hs_trend.xi_mu;
    xi_0 = input_dyn_hs_trend.xi_0;
    sig2_y = input_dyn_hs_trend.sig2_y;
    y = input_dyn_hs_trend.y;
    prior_phi = input_dyn_hs_trend.prior_phi; % Parameters for Beta prior
    H  = input_dyn_hs_trend.H;
    

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
        xi = 1./pgdraw(eta);
        iSig_xi = sparse(1:T,1:T,1./xi);

        % (c) Sample mu
        Kmu = xi_mu + xi_0 + (1-phi^2)*sum(xi(2:T));
        Cm = chol(Kmu,'lower');
        mu_hat = Kmu\((xi_mu + xi_0*h(1) + (1-phi)*sum(xi(2:T).*[h(2:T)-phi*h(1:T-1)])));
        mu = mu_hat + Cm'\randn(1,1);
        mu_tilde = [mu;(1-phi)*ones(T-1,1)*mu];

        % (d) Sample xi_mu
        xi_mu = 1./pgdraw(mu);

        % (e) Sample xi_0 
        xi_0 = 1./pgdraw(h_tilde(1));

        % (f) Sample phi (needs slice sampling)
            % define y_ar
            y_ar = h_tilde(2:end)./(xi(2:end).^(1/2));
            x_ar = h_tilde(1:end-1)./(xi(1:end-1).^(1/2));

            % define x_ar
        log_likelihood = @(x) -0.5*sum((y_ar - (2*x - 1).*x_ar).^2) + log(betapdf(x, prior_phi(1), prior_phi(2)));
        phi = slice_sampling(phi, log_likelihood, 0, 1); % Adjust upper and lower bounds as needed
        Hphi = speye(T) - sparse(2:T,1:(T-1),phi*ones(1,T-1),T,T);


    %% 2: Sample State variables
        % (a) tau
        iSig_tau = sparse(1:T,1:T,1./exp(h));
        Ktau= (eye(T)*1./sig2_y + H'*iSig_tau*H);
        Ctau = chol(Ktau,'lower');
        tau_hat = Ktau\(1./sig2_y*y);
        tau = tau_hat + Ctau'\randn(T,1);


    %% 3. Pack up samples
    output_dyn_hs_trend.tau = tau;
    output_dyn_hs_trend.h = h;
    output_dyn_hs_trend.phi = phi;
    output_dyn_hs_trend.h_tilde = h_tilde;
    output_dyn_hs_trend.Hphi = Hphi;
    output_dyn_hs_trend.iSig_xi = iSig_xi; 
    output_dyn_hs_trend.mu_tilde = mu_tilde;
    output_dyn_hs_trend.mu = mu;
    output_dyn_hs_trend.xi_mu = xi_mu;
    output_dyn_hs_trend.xi_0 = xi_0;
    output_dyn_hs_trend.sig2_y = sig2_y;
    output_dyn_hs_trend.y = y;
    output_dyn_hs_trend.prior_phi = prior_phi; % Parameters for Beta prior
    output_dyn_hs_trend.H = H;


end
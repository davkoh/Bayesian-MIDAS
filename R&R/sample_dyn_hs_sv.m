function [output_dyn_hs_sv] = sample_dyn_hs_sv(ystar,input_dyn_hs_sv)

% Define useful things
T = length(ystar);

% Unpack input
    g = input_dyn_hs_sv.g;
    phi = input_dyn_hs_sv.phi;
    g_tilde = input_dyn_hs_sv.g_tilde;
    Hphi = input_dyn_hs_sv.Hphi;
    iSig_xi = input_dyn_hs_sv.iSig_xi; 
    mu_tilde  = input_dyn_hs_sv.mu_tilde;
    mu = input_dyn_hs_sv.mu;
    xi_mu = input_dyn_hs_sv.xi_mu;
    xi_0 = input_dyn_hs_sv.xi_0;
    y = input_dyn_hs_sv.y;
    prior_phi = input_dyn_hs_sv.prior_phi; % Parameters for Beta prior
    H  = input_dyn_hs_sv.H;
    

% normal mixture
        pj = [.22715 .20674 .18842 .13057 .12047 .05591 .04775 .01575	.00609, .00115];
        mj = [-0.85173 .02266 -1.97278 .73504 -3.46788 -5.55246 1.34744 -8.68384 1.92677 -14.65];  % warning: means already adjusted
        sigj2 = [.62699 .40611 .98583 .26768 1.57469 2.54498 .17788 4.16591 .11265 7.33342];
        sigj = sqrt(sigj2);
        % sample S from a 10-point distrete distribution
        temprand = rand(T,1);
        q = repmat(pj,T,1).*normpdf(repmat(ystar,1,10),...
        repmat(g_tilde,1,10)+repmat(mj,T,1),repmat(sigj,T,1));
        q = q./repmat(sum(q,2),1,10);
        S = 10 - sum(repmat(temprand,1,10)<cumsum(q,2),2)+1;
        m_s = mj(S)'; iSig_v = sparse(1:T,1:T,1./sigj2(S));

        Kg= (iSig_v + Hphi'*iSig_xi*Hphi);
        Cg = chol(Kg,'lower');
        gtilde_hat = Kg\(iSig_v*(ystar- m_s -mu_tilde));
        g_tilde = gtilde_hat + Cg'\randn(T,1);
        g = g_tilde + mu;

        
        % (b) Sample xi_t (PG variables)
        eta = Hphi*g_tilde;
        xi = 1./pgdraw(eta);
        iSig_xi = sparse(1:T,1:T,1./xi);

        % (c) Sample mu
        Kmu = xi_mu + xi_0 + (1-phi^2)*sum(xi(2:T));
        Cm = chol(Kmu,'lower');
        mu_hat = Kmu\(( xi_0*g(1) + (1-phi)*sum(xi(2:T).*[g(2:T)-phi*g(1:T-1)])));
        mu = mu_hat + Cm'\randn(1,1);
        mu_tilde = [mu;(1-phi)*ones(T-1,1)*mu];

        % (d) Sample xi_mu
        xi_mu = 1./pgdraw(mu);

        % (e) Sample xi_0 
        xi_0 = 1./pgdraw(g_tilde(1));

        % (f) Sample phi (needs slice sampling)
            % define y_ar
            y_ar = g_tilde(2:end)./(xi(2:end).^(1/2));
            x_ar = g_tilde(1:end-1)./(xi(1:end-1).^(1/2));

            % define x_ar
            phi_temp = (phi + 1)/2 ;
            log_likelihood = @(x) -0.5*sum((y_ar - (2*x - 1).*x_ar).^2) + log(betapdf(x, prior_phi(1), prior_phi(2)));
            phi_temp= slice_sampling(phi_temp, log_likelihood, 0, 1); % Adjust upper and lower bounds as needed
            phi = 2*phi_temp -1;
            Hphi = speye(T) - sparse(2:T,1:(T-1),phi*ones(1,T-1),T,T);



    %% 3. Pack up samples
    output_dyn_hs_sv.g = g;
    output_dyn_hs_sv.phi = phi;
    output_dyn_hs_sv.g_tilde = g_tilde;
    output_dyn_hs_sv.Hphi = Hphi;
    output_dyn_hs_sv.iSig_xi = iSig_xi; 
    output_dyn_hs_sv.mu_tilde = mu_tilde;
    output_dyn_hs_sv.mu = mu;
    output_dyn_hs_sv.xi_mu = xi_mu;
    output_dyn_hs_sv.xi_0 = xi_0;
    output_dyn_hs_sv.y = y;
    output_dyn_hs_sv.prior_phi = prior_phi; % Parameters for Beta prior
    output_dyn_hs_sv.H = H;


end
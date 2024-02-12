% This function draws the log-volatilities in an unobserved components 
% in the noncentered parameterization

% Inputs: a0 and b0 are the prior mean and variance of h0;
%         Vomegah is the prior variance of omegah
function [h_tilde,h0,omegah,omegah_hat,Domegah,Vomegah] = ...
    SVRW_gam_hier(ystar,h_tilde,h0,omegah,a0,b0,Vomegah,V_a,V_b)

    T = length(h_tilde);
        % normal mixture
    pj = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
    mj = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819]...
        - 1.2704;  % warning: means already adjusted
    sigj2 = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
    sigj = sqrt(sigj2);
        % sample S from a 7-point distrete distribution
    temprand = rand(T,1);
    q = repmat(pj,T,1).*normpdf(repmat(ystar,1,7),...
        repmat(h0+omegah*h_tilde,1,7)+repmat(mj,T,1),repmat(sigj,T,1));
    q = q./repmat(sum(q,2),1,7);
    S = 7 - sum(repmat(temprand,1,7)<cumsum(q,2),2)+1;

        % sample h_tilde
    H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);    
    d_s = mj(S)'; iOs = sparse(1:T,1:T,1./sigj2(S));
    Kh = H'*H + omegah^2*iOs;
    h_tilde_hat = Kh\(iOs*omegah*(ystar-d_s-h0));
    h_tilde = h_tilde_hat + chol(Kh,'lower')'\randn(T,1);

        % sample h0 and omegah
    Xbeta = [ones(T,1) h_tilde];
    iVbeta = diag([1/b0 1/Vomegah]);    
    Kbeta = iVbeta + Xbeta'*iOs*Xbeta;
    beta_hat = Kbeta\(iVbeta*[a0;0] + Xbeta'*iOs*(ystar-d_s));
    beta = beta_hat + chol(Kbeta,'lower')'\randn(2,1);
    h0 = beta(1); omegah = beta(2);
        % randomly permute the signs h_tilde and omegah
    U = -1 + 2*(rand>0.5);
    h_tilde = U*h_tilde;
    omegah = U*omegah;

        % Sample from the posterior of Vomegah
    Vomegah = 1/gamrnd(V_a+1/2,1/(V_b+((omegah).^2)/2));

        % compute the mean and variance of the
        %    conditional density of omegah    
    Dbeta = Kbeta\speye(2);
    omegah_hat = beta_hat(2);
    Domegah = Dbeta(2,2);
end
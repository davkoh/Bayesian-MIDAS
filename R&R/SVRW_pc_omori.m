% This function draws the log-volatilities in an unobserved components 
% in the noncentered parameterization

% Inputs: a0 and b0 are the prior mean and variance of h0;
%         Vomegah is the prior variance of omegah
function [h_tilde,h0,omegah,omegah_hat,Domegah,Vomegah,Vh0] = ...
    SVRW_hier_pc(ystar,h_tilde,h0,omegah,a0,b0,Vomegah,Vh0)

    T = length(h_tilde);
        % normal mixture
    pj = [.22715 .20674 .18842 .13057 .12047 .05591 .04775 .01575	.00609, .00115];
    mj = [-0.85173 .02266 -1.97278 .73504 -3.46788 -5.55246 1.34744 -8.68384 1.92677 -14.65];  % warning: means already adjusted
    sigj2 = [.62699 .40611 .98583 .26768 1.57469 2.54498 .17788 4.16591 .11265 7.33342];
    sigj = sqrt(sigj2);
        % sample S from a 7-point distrete distribution
    temprand = rand(T,1);
    q = repmat(pj,T,1).*normpdf(repmat(ystar,1,10),...
        repmat(h0+omegah*h_tilde,1,10)+repmat(mj,T,1),repmat(sigj,T,1));
    q = q./repmat(sum(q,2),1,10);
    S = 10 - sum(repmat(temprand,1,10)<cumsum(q,2),2)+1;

        % sample h_tilde
    H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);    
    d_s = mj(S)'; iOs = sparse(1:T,1:T,1./sigj2(S));
    Kh = H'*H + omegah^2*iOs;
    h_tilde_hat = Kh\(iOs*omegah*(ystar-d_s-h0));
    
    try
        Ch = chol(Kh,'lower');
    catch
        Ch = chol(nearestSPD(full(Kh)),'lower'); %SVD does not support sparse matrix so convert Kbeta to normal first
        Ch = sparse(Ch); %Convert back to sparse matrix
    end


    h_tilde = h_tilde_hat + Ch'\randn(T,1);

        % sample h0 and omegah
    Xbeta = [ones(T,1) h_tilde];
    iVbeta = diag([1/Vh0 1/Vomegah]);    
    Kbeta = iVbeta + Xbeta'*iOs*Xbeta;
    beta_hat = Kbeta\(iVbeta*[a0;0] + Xbeta'*iOs*(ystar-d_s));

    try
        Cbeta = chol(Kbeta,'lower');
    catch
        Cbeta = chol(nearestSPD(full(Kbeta)),'lower'); %SVD does not support sparse matrix so convert Kbeta to normal first
        Cbeta= sparse(Cbeta); %Convert back to sparse matrix
    end



    beta = beta_hat + Cbeta'\randn(2,1);
    h0 = beta(1); omegah = beta(2);
        % randomly permute the signs h_tilde and omegah
    U = -1 + 2*(rand>0.5);
    h_tilde = U*h_tilde;
    omegah = U*omegah;

        % Sample from the posterior of Vomegah
    Vomegah = sample_V2_slice(Vomegah,omegah,0,0,0.5,60);
        % Sample from the posterior of Vomegah
    Vh0 = sample_V2_slice(Vh0,h0,0,0,0.5,60);

        % compute the mean and variance of the
        %    conditional density of omegah    
    Dbeta = Kbeta\speye(2);
    omegah_hat = beta_hat(2);
    Domegah = Dbeta(2,2);
end
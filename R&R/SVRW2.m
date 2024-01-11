
function h = SVRW2(Ystar,h,sig,h0)

T = length(h);
    % define normal mixture
pi = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mi = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;  %% means already adjusted!! %%
sigj2 = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
sqrtsigi = sqrt(sigj2);

    % sample S from a 7-point distrete distribution
tmprand = rand(T,1);
q = repmat(pi,T,1).*normpdf(repmat(Ystar,1,7),repmat(h,1,7)+repmat(mi,T,1), repmat(sqrtsigi,T,1));
q = q./repmat(sum(q,2),1,7);
S = 7 - sum(repmat(tmprand,1,7)<cumsum(q,2),2)+1;
    
    % sample h
Hh = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
iSh = speye(T).*1./sig;
dconst = mi(S)'; iOmega = sparse(1:T,1:T,1./sigj2(S));
%alph = Hh\[h0;sparse(T-1,1)]; 
HiSH_h = Hh'*iSh*Hh;
Kh = HiSH_h + iOmega;
h_hat = Kh\(HiSH_h*h0*ones(T,1)+ iOmega*(Ystar-dconst));
h = h_hat + chol(Kh,'lower')'\randn(T,1);
end
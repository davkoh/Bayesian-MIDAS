function [h,phih,Sigmah,muh,lambdah,tauh] = sample_sv_hs_ar1(s2,h,sv_prior,phih,Sigmah,muh,lambdah,tauh)
    %SAMPLE_SV_AR1 Summary of this function goes here
%   Detailed explanation goes here

% Unpack the SV_prior
phih0 = sv_prior.phih0;
Vphih = sv_prior.Vphih;
muh0 = sv_prior.muh0;
Vmuh = sv_prior.Vmuh;

%% normal mixture
pimix = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mi = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;  %% means already adjusted!! %%
sigi = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
sqrtsigi = sqrt(sigi);
T = length(h);


%% sample h

Ystar = log(s2 + 0.001);

Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);

HiSH = Hphi'*sparse(1:T,1:T,[(1-phih^2)/Sigmah(1); 1./Sigmah(2:end).*ones(T-1,1)])*Hphi;

deltah = Hphi\[muh; muh*(1-phih)*ones(T-1,1)];

%% sample S from a 7-point distrete distribution
temprand = rand(T,1);
qns = repmat(pimix,T,1).*normpdf(repmat(Ystar,1,7),repmat(h,1,7)+repmat(mi,T,1), repmat(sqrtsigi,T,1));
q = qns./repmat(sum(qns,2),1,7);
S = 7 - sum(repmat(temprand,1,7)<cumsum(q,2),2)+1;

%% sample h
muLarge = mi(S)'; isiglarge = 1./sigi(S);
invSigS = sparse(1:T,1:T,isiglarge);
Kh =  HiSH + invSigS;
hhat = Kh\(HiSH*deltah+invSigS*(Ystar-muLarge));
h = hhat + chol(Kh,'lower')'\randn(T,1);

%% sample Sigmah
errh = [(h(1)-muh)*sqrt(1-phih^2);  h(2:end)-phih*h(1:end-1)-muh*(1-phih)];
v=1./gamrnd(1,1./(1+1./lambdah));
lambdah=1./gamrnd(1, 1./(1./v + (errh.^2)./(2.*tauh)));
xi=1./gamrnd(1,1./(1+1./tauh));
tauh=1./gamrnd((T+1)/2, 1./(1./xi +0.5*sum(sum(errh.^2./lambdah))  ));
Sigmah=(lambdah.*tauh)+1e-10;


%% sample phih
iS = 1./(Sigmah.*ones(T,1));
Xphi = h(1:end-1)-muh;
yphi = h(2:end) - muh;
Dphi = 1/(1/Vphih + (Xphi.*iS(2:end))'*Xphi);
phihat = Dphi*(phih0/Vphih + (Xphi.*iS(2:end))'*yphi);
phic = phihat + sqrt(Dphi)*randn;
g = @(x) -.5*log((1./iS(1))./(1-x.^2))-.5*(1-x.^2)*iS(1)*(h(1)-muh)^2;
if abs(phic)<.9999
    alpMH = exp(g(phic)-g(phih));
    if alpMH>rand
        phih = phic;
    end
end
%% sample muh
Xm = ones(T-1,1)*(1-phih);
Dmuh = 1/(1/Vmuh + ( (Xm.*iS(2:end))'*Xm + (1-phih^2).*iS(1)) ); 
Ym = h(2:end)-phih*h(1:end-1);
muhhat = Dmuh*(muh0/Vmuh + (1-phih^2)*iS(1)*h(1) + (Xm.*iS(2:end))'*Ym);
muh = muhhat + sqrt(Dmuh)*randn;
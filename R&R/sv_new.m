%% Updating the SV Sampler %%
clear all;

% Choose Parameters of the DGP
T = 200;
sd_tau = 0.1;
t_0 = 0;
sd_y = 1;
h_0 = 0;
sd_h = 0.5;
phih = 0.95;

% Generate Data
tau = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);
outlier = 0;

for t= 1:T
    if t == 1
    
        h(1) = h_0 + sd_h*randn;
    else
    
    h(t) = phih*h(t-1) + sd_h*randn;
    end

    if t == 70
    
        h(t) = h(t);
        y(t) =  exp(1/2*h(t))*randn -outlier;
    elseif t ==71
        h(t) = h(t);
        y(t) =  exp(1/2*h(t))*randn +outlier;
    else
        y(t) = exp(1/2*h(t))*randn;
    end
end

h_true = h;

clf;
plot(y)
plot(h_true)


% Priors
V_t0 = 10; % Variance of initial condition
nu_sig0 = 3; S_sig0 = 1*(nu_sig0-1); % \sigma^2
Vh = 10;
nu_omega0 = 3; S_omega0 = .25^2*(nu_omega0-1); % \omega_tau^2

% Set up difference matrix
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);

% Initalise series
lam = ones(T,1);
nu = 1;
eta_lam = ones(T,1);
eta_nu = 1;
tau0 = 0;
Sigma_tau = ones(T,1);
sig2 = 1;
tauh=1;
lambdah=ones(T,1);
Sigh=(lambdah.*tauh);
h = ones(T,1);
h0 = log(var(y));
omegah2 = 0.1;
muh = 0;

 omega_tau2 = .1;

% Sampler Specs
rng(1,'twister');
iter = 20000;
burn = iter/2;
mcmc = iter-burn;

%% Gibbs Sampler

% Storage
tau_save = NaN(T,mcmc);
tau0_save = NaN(mcmc,1);
sig2_save = NaN(mcmc,1);
lam_save = NaN(T,mcmc);
nu_save = NaN(mcmc,1);
h_save = NaN(T,mcmc);
h0_save = NaN(mcmc,1);

mu_save = NaN(mcmc,1);
phih_save = NaN(mcmc,1);
omegah_save = NaN(mcmc,1);



pr_S_sigh = 0.1; % Prior on the state variation;
pr_nu_sigh = 5;

phih0 = .97; Vphih = .1^2;

muh0 = 1; Vmuh = 10;




for i = 1:iter

s2 = (y).^2;


    %SAMPLE_SV_AR1 Summary of this function goes here
%   Detailed explanation goes here
%% normal mixture
pimix = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mi = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;  %% means already adjusted!! %%
sigi = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
sqrtsigi = sqrt(sigi);
T = length(h);


%% sample h

Ystar = log(s2 + 0.001);

Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);

HiSH = Hphi'*sparse(1:T,1:T,[(1-phih^2)/omegah2; 1/omegah2*ones(T-1,1)])*Hphi;

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

%% sample omegah2
errh = [(h(1)-muh)*sqrt(1-phih^2);  h(2:end)-phih*h(1:end-1)-muh*(1-phih)];
omegah2 = 1/gamrnd((length(errh) + pr_nu_sigh)/2, 2/(sum(errh.^2) + pr_S_sigh));


%% sample phih
iS = 1./(omegah2*ones(T,1));
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

    if i> burn 
    isave = i - burn;
    lam_save(:,isave) = lam;
    omegah_save(isave) = sqrt(omegah2);
    h_save(:,isave) = h;
    mu_save(isave) = muh;
    phih_save(isave) = phih;
end


end



% Evaluation
ph = mean(h_save,2);
clf;
plot(exp(0.5*ph))
hold on
plot(exp(0.5*h_true))

hhat = mean(exp(h_save/2),2);  %% plot std dev
plot(hhat)









%% Copied over function from the dfm code

function [h,lamh,nu,muh,phih,omegah2] = sv_ar1(s2,h,lamh,nu,muh,phih,omegah2,phih0,Vphih,muh0,Vmuh,pr_nu_sigh,pr_S_sigh,nu_ub)
%SAMPLE_SV_AR1 Summary of this function goes here
%   Detailed explanation goes here
%% normal mixture
pimix = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mi = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;  %% means already adjusted!! %%
sigi = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
sqrtsigi = sqrt(sigi);
T = length(h);


%% sample h

Ystar = log(s2 + 0.001);

Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);

HiSH = Hphi'*sparse(1:T,1:T,[(1-phih^2)/omegah2; 1/omegah2*ones(T-1,1)])*Hphi;

deltah = Hphi\[muh; muh*(1-phih)*ones(T-1,1)];

%% sample S from a 7-point distrete distribution
temprand = rand(T,1);
qns = repmat(pimix,T,1).*normpdf(repmat(Ystarnonan,1,7),repmat(h,1,7)+repmat(mi,T,1), repmat(sqrtsigi,T,1));
q = qns./repmat(sum(qns,2),1,7);
S = 7 - sum(repmat(temprand,1,7)<cumsum(q,2),2)+1;

%% sample h
muLarge = mi(S)'; isiglarge = 1./sigi(S);
invSigS = sparse(1:T,1:T,isiglarge);
Kh =  HiSH + invSigS;
hhat = Kh\(HiSH*deltah+invSigS*(Ystar-muLarge));
h = hhat + chol(Kh,'lower')'\randn(T,1);

%% sample omegah2
errh = [(h(1)-muh)*sqrt(1-phih^2);  h(2:end)-phih*h(1:end-1)-muh*(1-phih)];
omegah2 = 1/gamrnd((length(errh) + pr_nu_sigh)/2, 2/(sum(errh.^2) + pr_S_sigh));


%% sample phih
iS = 1./(omegah2*ones(T,1));
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


    Ss = nu*ones(T,1);
    Ss(idx_nan) = Ss+s2.*exp(-h);
    nus = nu*ones(T,1);
    nus(idx_nan) = nus(idx_nan) + 1;
    lamh = 1./gamrnd(nus./2,2./Ss);
    nu = sample_nu(lamh,nu,nu_ub);

end


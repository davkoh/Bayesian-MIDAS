% Bayesian Adaptive Group Lasso with Spike-and-Slab

function out=BMIDAS_AGLasso_SS_independent(x,y,X_fore,It,groups)

warning('off','all')

T = size(x,1);
p = size(x,2);
ng = size(unique(groups),2); % number of groups;
nj = histc(groups, unique(groups));
%nj = ones(ng,1)*size(Q,1);   % number of predictors in each group
groups = cell(ng,1);
groups{1}=[1:nj(1)];
for i=2:ng
    groups{i}=[sum(nj(1:i-1))+1:sum(nj(1:i-1))+nj(i)];
end

x1 = x;
[x1,mu_x,sig_x] = standardize(x);
mu_x(:,:) = 0;
sig_x(:,:) = 1;

[y,mu_y,sig_y] = center(y);

X_fore_star=(X_fore-mu_x)./sig_x;

xg = cell(ng,1);
sg = NaN(ng,1);
for ig = 1:ng
    xg{ig} = x(:,groups{ig});
    sg(ig) = length(groups{ig});
end
  
% Gibbs parameters
nsave = It.nsave;
nburn = It.nburn;
nthin = It.nthin;
ntot = nsave + nburn;
ndraw = nsave/nthin;

beta_draws = zeros(ndraw,p+1);
tau2_draws = zeros(ndraw,ng);
lambda2_draws = zeros(ndraw,ng);
sig2_draws = zeros(ndraw,1);
y_fore_draws = zeros(ndraw,1);
pi0_draws = zeros(ndraw,1);
pi1_draws = zeros(ndraw,ng);
Z_draws = zeros(ndraw,ng);

% initial values for stabilization algorithm
zeta=1*ones(1,ng);         
kappa=0;
nu=0;

logmax = log(3);
logmin = log(1);
diff_ = (logmin-logmax)/(ntot-1);
range = exp(logmax:diff_:logmin);
eps_zeta = exp(logmax)*ones(1,ng);

% Initialization
XgtXg = cell(ng,1);
for ig=1:1:ng
    XgtXg{ig} = x(:,groups{ig})'*x(:,groups{ig});
end

% Priors
lambda2 = ones(1,ng);
sig2 = 1;
tau2 = NaN(1,ng);
beta = cell(ng,1);
betan = NaN(ng,1);
for ig=1:ng
    tau2(1,ig) = (1./(lambda2(1,ig)/2)).*randg((sg(ig)+1)./2,[1 1]) + 1e-10;
    cov = eye(sg(ig))*(sig2*tau2(1,ig));
    beta{ig} = zeros(1,sg(ig)) + randn(1,sg(ig))*chol(cov);
end
beta=(cell2mat(beta'))';

kappabar=(1+1/ng);
u=kappabar;
aa = kappabar*(ng^u);
bb = 1;
g1 = randg(aa,[1 1]);
g2 = randg(bb,[1 1]);
pi0 = g1./(g1 + g2);
pi1 = zeros(1,ng);
Z = zeros(1,ng);
nu0 = 3; S0 = 1*(nu0 - 1);

%% GIBBS ITERATIONS

for irep = 1:ntot
           
    % 1. Update beta and tau2
    u_rand = rand(ng,1);
    for ig = 1:ng
        D = (1/tau2(ig)).*eye(sg(ig)); 
        A = XgtXg{ig}/sig2+D; 
        AA = (A + A.')/2;
        beta(groups{ig})=zeros(sg(ig),1);
        Z(ig) = 0;
        betan(ig) = 0;
        b = y - x*beta;
        xb = 1/sig2*x(:,groups{ig})'*b;

        maxAA = max(AA(:));
        L = (-sg(ig)/2)*log(tau2(ig)) + (-1/2)*log(det(AA/maxAA)) + (-size(AA,1)/2)*log(maxAA) + (0.5/sig2)*(xb'*(AA\xb));
        pi1(ig) = pi0/(pi0+(1-pi0)*exp(L));

        if u_rand(ig)>=pi1(ig)
           beta(groups{ig}) = ((AA\xb)' + randn(1,sg(ig))*chol((AA\eye(sg(ig)))))';
           Z(ig) = 1;
           betan(ig) = norm(beta(groups{ig}),2);
           a1 = sqrt(lambda2(ig))./betan(ig);
           a2 = lambda2(ig);
           tau2(ig) = (1./randig(a1,a2)) + 1e-10; 
        end       
    end
    tau2(Z==0) = (1./(lambda2(Z==0)./2)).*randg((sg(Z==0)'+1)./2,[1 size(lambda2(Z==0),2)]) + 1e-10;
    
    % 2. Update sigma2 from Inverse Gamma
    sig2 = 1/gamrnd(nu0+T/2,1/(S0 + (y-x*beta)'*(y-x*beta)/2));

    % 3. Update pi from Beta
    g1 = randg(aa+ng-sum(Z),[1 1]);
    g2 = randg(bb+sum(Z),[1 1]);
    pi0 = g1./(g1 + g2); 
    
    % 4. Update lambda2_j
    a_n=1./(zeta.^0.8);
    s_i_1=log(sqrt(lambda2));
    s_i=s_i_1+a_n.*((sg'+1)-exp(2*s_i_1).*tau2);
    lambda2=exp(2*s_i);
    % Stabilization algorithm        
    Kappa=ones(ng,1)*[max([-5 (-kappa-1)]),kappa+1];
    Delta=abs(s_i'-s_i_1')';
    if mean(s_i'>=Kappa(:,1))==1 && mean(s_i'<=Kappa(:,2))==1 && mean(Delta<=eps_zeta)==1
       zeta(Z==1)=zeta(Z==1)+1;
       nu=nu+1;
       eps_zeta(Z==1) = exp(log(eps_zeta(Z==1))+diff_);
    else
       ii=1;
       while (mean(s_i'<Kappa(:,1))>0 || mean(s_i'>Kappa(:,2))>0) || mean(Delta>eps_zeta)>0
             l1=min([s_i_1(s_i>=Kappa(:,2)'); Kappa(s_i'>=Kappa(:,2),2)']);
             u1=max([s_i_1(s_i>=Kappa(:,2)'); Kappa(s_i'>=Kappa(:,2),2)']);
             s_i(1,s_i'>=Kappa(:,2))=l1 + (u1-l1).*rand(1,length(l1));
             l2=min([s_i_1(s_i<=Kappa(:,1)'); Kappa(s_i'<=Kappa(:,1),1)']);
             u2=max([s_i_1(s_i<=Kappa(:,1)'); Kappa(s_i'<=Kappa(:,1),1)']);             
             s_i(1,s_i'<=Kappa(:,1))=l2 + (u2-l2).*rand(1,length(l2));
             lambda2=exp(2*s_i);
             Delta=abs(s_i'-s_i_1')';
             tau2 = (1./(lambda2./2)).*randg((sg'+1)./2,[1 size(lambda2,2)]) + 1e-10;
             ii=ii+1;
             if ii>=100
                break
             end
       end
       zeta(Z==1)=zeta(Z==1)+1;
       kappa=kappa+1;
       nu=0;
       eps_zeta(Z==1) = exp(log(eps_zeta(Z==1))+diff_);
    end
    
    if ((irep>nburn) && (mod(irep - nburn, nthin) == 0))
       intercept = mu_y + sqrt(sig2/T)*randn;
       beta_draws((irep-nburn)/nthin,:) = [intercept beta'];
       tau2_draws((irep-nburn)/nthin,:) = tau2;
       sig2_draws((irep-nburn)/nthin) = sig2; 
       lambda2_draws((irep-nburn)/nthin,:) = lambda2;
       pi0_draws((irep-nburn)/nthin) = pi0;
       Z_draws((irep-nburn)/nthin,:) = Z;
       pi1_draws((irep-nburn)/nthin,:) = pi1;
       y_fore_draws((irep-nburn)/nthin,:) = ([1 X_fore_star]*[intercept;beta]) + sqrt(sig2)*randn;
    end   

end


%% Output

out.beta0 = beta_draws(:,1);
temp = arrayfun(@(k) ((beta_draws(:,groups{k}+1)./sig_x(1,groups{k}))')',1:ng,'UniformOutput',false)';
%temp = arrayfun(@(k) ((beta_draws(:,groups{k}+1)./sig_x(1,groups{k}))')',1:ng,'UniformOutput',false)';
out.beta = cat(2,temp{:});
out.lambda2 = lambda2_draws;
out.tau2 = tau2_draws;
out.sigma2 = sig2_draws;
out.Z = median(Z_draws);
out.y_fore = y_fore_draws;
out.y_fore_median = median(y_fore_draws);    % median estimator for SS




end

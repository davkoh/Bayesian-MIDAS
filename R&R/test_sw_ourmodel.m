%% Generate from SW 2007 model to see whether the HS model can perform adaquately %%

clear all;
rng(1,'twister');

% Generate Data


% Choose Parameters of the DGP
T = 200;
sd_tau = 0.1;
t_0 = 0;
sd_y = 1;
h_0 = 0;
sd_h = 0.1;
g_0 = 0;
sd_g = 0.1;

% Generate Data
tau = zeros(T,1);
g = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);

outlier = 150;


for t= 1:T
    if t == 1
        tau(1) = t_0 + sd_tau*randn;
        h(1) = h_0 + sd_h*randn;
        g(1) = g_0 + sd_g*randn;
    else
    g(t) = g(t-1) + sd_g*randn;
    tau(t) = tau(t-1) + exp(1/2*g(t))*randn;
    h(t) = h(t-1) + sd_h*randn;
    end

    
    y(t) = tau(t) + exp(1/2*h(t))*randn;
    
end

tau_true = tau;
h_true = h;
g_true = g;

clf;
plot(y)
plot(h_true)
plot(tau_true)
plot(g_true)


nsim = 10000;
burnin = 10000;
valh = 0;
valg = 0;

%% prior
a0_h = 0; b0_h = 10;
a0_g = 0; b0_g = 10;
a0_tau = 0; b0_tau = 10;
Vomegah = .0001;
Vomegag = .0001;
    
% initialize the Markov chain
h0 = log(var(y))/5; g0 = log(var(y))/10; tau0 = mean(y);
omegah = sqrt(.2);
omegag = sqrt(.2);
h_tilde = zeros(T,1);
g_tilde = zeros(T,1);
h = h0 + omegah*h_tilde;
g = g0 + omegah*g_tilde; 

% define a few things
n_grid = 500; % number of grid points
omh_grid = linspace(-1,1,n_grid)';
omg_grid = linspace(-1,1,n_grid)';
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);

% initialize for storeage
store_theta = zeros(nsim,5); % [omegah omegag h0 g0 tau0]
store_tau = zeros(nsim,T); 
store_h = zeros(nsim,T);
store_g = zeros(nsim,T);
store_pomh = zeros(n_grid,1);
store_pomg = zeros(n_grid,1);
store_lpden = zeros(nsim,3); % log posterior densities of 
                             % omegah = valh; omegag = valg;
                             % omegah = omegag 
                            
rand('state', sum(100*clock) ); randn('state', sum(200*clock) );    
for isim = 1:nsim+burnin
    
        % sample tau    
    iOh = sparse(1:T,1:T,1./(exp(h).*lamt));
    HiOgH = H'*sparse(1:T,1:T,1./exp(g))*H;
    Ktau =  HiOgH + iOh;    
    tau_hat = Ktau\(tau0*HiOgH*ones(T,1) + iOh*y);
    tau = tau_hat + chol(Ktau,'lower')'\randn(T,1);
    
        % sample h_tilde 
    ystar = log((y-tau).^2.*lamt.^(-1) + .0001);
    [h_tilde h0 omegah omegah_hat Domegah] = ...
        SVRW_gam(ystar,h_tilde,h0,omegah,a0_h,b0_h,Vomegah); 
    h = h0 + omegah*h_tilde;    
    
        % sample g_tilde
    ystar = log((tau-[tau0;tau(1:end-1)]).^2 + .0001);
    [g_tilde g0 omegag omegag_hat Domegag] = ...
        SVRW_gam(ystar,g_tilde,g0,omegag,a0_g,b0_g,Vomegag); 
    g = g0 + omegag*g_tilde;
    
        % sample tau0
    Ktau0 = 1/b0_tau + 1/exp(g(1));
    tau0_hat = Ktau0\(a0_tau/b0_tau + tau(1)/exp(g(1)));
    tau0 = tau0_hat + sqrt(Ktau0)'\randn;

        %% Sample t-errors
    [lamt,nut] = sample_t_errors(y-tau,lamt,nut,nub,exp(h));
            
    if (mod(isim, 1000) == 0)
        disp([num2str(isim) ' loops... ']);
    end     
    
    if isim > burnin
        isave = isim - burnin;
        store_tau(isave,:) = tau';
        store_h(isave,:) = h'; 
        store_g(isave,:) = g'; 
        store_theta(isave,:) = [omegah omegag h0 g0 tau0]; 
        save_nut_sw(isave) = nut;
            score_save(isave) = mean(normpdf(y,tau,lamt.*exp(h)));


        lh0 = -.5*log(2*pi*Domegah) - .5*(omegah_hat-valh)^2/Domegah;
        lg0 = -.5*log(2*pi*Domegag) - .5*(omegag_hat-valg)^2/Domegag;
        lhg0 = lh0 + lg0;
        store_lpden(isave,:) = [lh0 lg0 lhg0];
        
        store_pomh = store_pomh + normpdf(omh_grid,omegah_hat,sqrt(Domegah));
        store_pomg = store_pomg + normpdf(omg_grid,omegag_hat,sqrt(Domegag));        
    end    
end

theta_hat = mean(store_theta)';
tau_hat = mean(store_tau)';
h_hat = mean(exp(store_h/2))'; 
g_hat = mean(exp(store_g/2))'; 

clf;
plot(tau_hat)
hold on
plot(tau_true)

clf;
plot(exp(0.5*h_hat))
hold on
plot(exp(0.5*h_true))
plot(exp(0.5*g_true))
plot(exp(0.5*h_true) + exp(0.5*g_true))

trendrmse = sqrt(sum((tau_hat - tau_true).^2)/T)
svrmse = sqrt(sum((h_hat - h_true).^2)/T)



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
sd_h = 0.01;
g_0 = 0;
sd_g = 0.01;

% Generate Data
tau = zeros(T,1);
g = zeros(T,1);
y = zeros(T,1);
h = zeros(T,1);

outlier = 20;


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

    if t == 70
        y(t) = tau(t) + exp(1/2*h(t))*randn - outlier;
    elseif t == 71
        y(t) = tau(t) + exp(1/2*h(t))*randn + outlier;
    else
        y(t) = tau(t) + exp(1/2*h(t))*randn;
    end
    
end

tau_true = tau;
h_true = h;
g_true = g;

clf;
plot(y)
plot(h_true)
plot(tau_true)
plot(g_true)

% Priors
V_t0 = 10; % Variance of initial condition
Vh = 10;

prior_tau.V_t0 = 10;
prior_tau.H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);


prior_sv.Vh = 10;



% Initalise series
lambdatau = ones(T,1);
nutau = 1;
eta_lam = ones(T,1);
eta_nu = 1;
tau0 = 0;
Sigmatau = ones(T,1);
sig2 = 1;
nuh=1;
lambdah=ones(T,1);
Sigmah=(lambdah.*nuh);
h = ones(T,1);
h0 = log(var(y));
tau= zeros(T,1);
nub = 50;
lamt = ones(T,1);
nut = 5;

% Sampler Specs
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
nut_save = NaN(mcmc,1);

for i = 1:iter
    %% Sample tau
    [tau,Sigmatau,tau0,lambdatau,nutau] = sample_trend_unitroot_hs(y,tau,Sigmatau,tau0,lambdatau,nutau,exp(h).*lamt,prior_tau);

    %% Sample h
    u = y-tau;
    Ystar = log(u.^2.*lamt.^(-1) + .0001);
    [h,Sigmah,h0,lambdah,nuh] = sample_sv_unitroot_hs(Ystar,h,Sigmah,h0,lambdah,nuh,prior_sv);

    %% Sample t-errors
    [lamt,nut] = sample_t_errors(u,lamt,nut,nub,exp(h));
    

if i> burn 
    isave = i - burn;
    tau_save(:,isave) = tau;
    tau0_save(isave) = tau0;
    lam_save(:,isave) = lambdatau;
    nu_save(isave) = nutau;
    h_save(:,isave) = h;
    h0_save(isave) = h0;
    nut_save(isave) = nut;
    score_save(isave) = mean(normpdf(y,tau,lamt.*exp(h)));
end

end

log_score = mean(score_save)

ptau = mean(tau_save,2);
clf;
plot(ptau)
hold on
plot(tau_true)

ph = mean(h_save,2);
clf;
plot(exp(0.5*ph))
hold on
plot(exp(0.5*h_true))
plot(exp(0.5*g_true))
plot(exp(0.5*h_true) + exp(0.5*g_true))

trendrmse = sqrt(sum((ptau - tau_true).^2)/T)
svrmse = sqrt(sum((ph - h_true).^2)/T)

%% Plots to show in the overleaf
clf; 
fig1 = figure;
plot(y,LineWidth=2)
title("Generated Target")
xlabel("Time")

% Comparing Trends
fig2 = figure;
plot(tau_true,LineWidth=2)
title("Comparing Trends")
hold on
plot(ptau,LineWidth=2)
plot(tau_hat,LineWidth=2)
xlabel("Time")
legend("Target","HS (RMSE = 0.49)","Our Model (RMSE = 0.56)")

% Comparing SVs
fig3 = figure;
plot(h_true,LineWidth=2)
title("Comparing SVs")
hold on
plot(ph,LineWidth=2)
plot(h_hat,LineWidth=2)
xlabel("Time")
legend("Target","HS (RMSE = 0.6)","Our Model without t-errors (RMSE = 1.58)")

% Comparing degrees of freedom
fig4 = figure;
hist(nut_save)
hold on
histogram(save_nut_sw,"FaceColor","black")
xline(50,LineWidth=2,Color="red")
legend("HS","Our Model","True Value (dof = \infty)")

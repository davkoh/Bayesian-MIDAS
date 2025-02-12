%% Estimation script for "Flexible Bayesian MIDAS: time‑variation, group‑shrinkage and sparsity"
%%% David Kohns, Aalto University
%%% Galina Potjagailo, Bank of England
%%% Date 12.02.2025

%%% This code estimates the T-SVt-BMIDAS model with GIGG prior and ex-post
%%% group sparsification
%   - This code has been optimised for the cluster and is currently under
%   construction to be more easily modifiable without
%%%% ------------------ %%%
clear all
rng(1,'twister');  %set seed



%% Define directory and upload data

mkdir 'Output'
outputfolder = char([cd,'\Output']);  
addpath("../Data")
addpath("../Matlab")

dat_choice =  dat_choice_bash; % Whether to use old data or new data

if dat_choice == 1
load("UK_dat_2024.mat");
else
load("matlab_input.mat");
end

calendar_mod

% Output Matrices
crps_all = zeros(vint,nfor); % Storage for CRPS values
y_pred_all = zeros(MCMC,vint,nfor); % stores predictive distributions for each nowcast
rtresid_all = zeros(vint,nfor); % stores residuals for each nowcast, based on mean predictive
rtlogscores_all = zeros(vint,nfor); % stores log-scores for each nowcast
yf_all =zeros(nfor,1); % Saves the out-of-sample LHS variable for each quarter
dq_nfor = []; % saves dates of quarters to be nowcasted
pincl = zeros(vint,max(unique(groupall)),nfor);
modall = zeros(vint,MCMC,nfor);

%% Define things for cluster
initParPool()

hyperpars = [1/T,1/T;1/T,0.5;1/T,1;1,1/T;0.5,1/T;1,1;0.5,0.5]; % Hyperpriors for the GIGG prior

trend = 1; % trend_bash
SV = 1; % sv_bash
t = 0; % t_bash
almonrest = 1; % almonrest_bash
group_sparse = 1; % group_sparse_bash
ortho_choice = 1; % ortho_choice_bash

for gg = 2:2 % Loop over all hyperparameters for the GIGG prior
tic
%% Loop over time periods
parfor (tperiod = 1:nfor, str2double(getenv('SLURM_CPUS_PER_TASK')))
Xf = (Xm(1:tin+tperiod,:));
yf = y(tin+tperiod:end,:); 
T=tin-1+tperiod;

    predv = []; % Locals for the scores
    rtscores = [];
    crpsv = [];
    
    pincl_temp = zeros(vint,max(unique(groupall))); 
    mod_temp = zeros(vint,MCMC);
            
yf_all(tperiod,:) = yf(1,1);

%% Loop over nowcast periods
for v = 1:vint
    if pseudo_cal == 1
       display (['Draws for hyperparameter ' num2str(gg) ', evaluation for time period ' num2str(tperiod) ' of ' num2str(nfor) ', nowcast period ' num2str(v) ' of ' num2str(vint)])    %%% one additional dimension here over which draws are looped?
    else
        display (['Draws for hyperparameter ' num2str(gg), ', evaluation for time period ' num2str(tperiod) ' of ' num2str(nfor)])    %%% one additional dimension here over which draws are looped?
    end
        
% Find lags
xind =  find(puball(v,:)==1); % Find index of lags to include depending on the publication calendar

% Create GIGG readable groups from UMIDAS data
grp_idx = NaN(size(xind,2),1); % Create artificial group index from 1:# of groups available (gigg code assumes equidistant group numbering)
num_grp = histc(groupall(xind), unique(groupall(xind)));
sum_grp = size(unique(groupall(xind)),2);
for ii = 1:sum_grp
    end_idx = sum(num_grp(1:ii));
    start_idx = end_idx - num_grp(ii) +1;
    grp_idx(start_idx:end_idx) = repmat(ii,num_grp(ii),1) ;
end

% Hardcoded adjustment for direct forecast of MIDAS component
if v <7
Xv = (Xm(1:tin-4+tperiod,xind));
T=tin-4+tperiod;
else
    Xv = (Xm(1:tin-1+tperiod,xind));
    T=tin-1+tperiod;
end

if almonrest == 1
% Transform UMIDAS to Almon MIDAS 
[Xv,Xall] = midas_dat_r2_final(Xv,grp_idx,poly); % 
Xv = (Xv);

% Get grouping after Almon Transformation
lagl = histc(grp_idx, unique(grp_idx)); % lag length of each covariate;
   grp_idx_temp = [];
   
   for ii = 1:size(lagl,1)
       if lagl(ii)>2
       tt = repmat(ii,2,1)';
       else
        tt=    repmat(ii,lagl(ii),1)';
       end
       grp_idx_temp = [grp_idx_temp tt];
   end
else 
    grp_idx_temp = grp_idx';
end

% Orthonormalise the data (QR decomposition) for validity of
% Group-sparsification in steps below
if ortho_choice == 1
%
Xsvd = NaN(T,size(Xv,2));
Qj = cell(sum_grp,1);
Lam_inv_sqr = cell(sum_grp,1);
for j= 1:sum_grp
    xind1 = find(grp_idx_temp == j);
    gramm_m = 1/T * Xv(:,xind1)'*Xv(:,xind1);
    [U,S,D] = svd(gramm_m);
    Qj{j} = U;
    Lam_inv_sqr{j} = diag(diag(S).^(-1/2));
    Xsvd(:,xind1) = Xv(:,xind1)*Qj{j}*Lam_inv_sqr{j}/sqrt(T);
end
Xv = Xsvd;
end


%%%%%%%%%%%%%%%%%%  Estimate Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input = [];
input.grp_idx = grp_idx_temp';
if v < 7
    input.Y = y(3:tin-2+tperiod,:);
else
    input.Y = y(1:tin-1+tperiod,:);
end
input.X = Xv;
input.burnin = BURNIN;
input.samples = MCMC;
input.btrick = 0;
input.a = repmat(hyperpars(gg,1),sum_grp,1);
input.b = repmat(hyperpars(gg,2),sum_grp,1);
input.standardise = 0;
input.trend_ind = trend;
input.sv_ind = SV;
input.t_ind = t ;


[out] = gigg_hier_bsts_sv_t_aut(input);


% Perform group sparsification
if group_sparse == 1
[beta_out] = group_savs_orth(out.beta,grp_idx_temp');
else
    beta_out= out.beta;
    betas_final = out.beta;
end

%  Transform back to non-orthogonalised
if ortho_choice == 1
betas_final = out.beta;
for j = 1:sum_grp
xind1 = find(grp_idx_temp == j);
betas_final(xind1,:) = Qj{j}*Lam_inv_sqr{j}*beta_out(xind1,:)/sqrt(tin);
end
end

% Variable Selection Info
idx_first_memb = [];
for iii = 1:size(unique(grp_idx_temp),2)
    idx_first_memb = [idx_first_memb;min(find(iii== grp_idx_temp))];
end
pincl_temp(v,unique(groupall(xind))) = (sum(betas_final(idx_first_memb,:)'~=0)/MCMC)' ;


%%%%%%%%%%%%  Nowcasting %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Retrieve parameters
betas_final = betas_final';
gout = out.g;
hout = out.h;
tauout = out.tau;
thetaout = out.theta;
nuout = out.nu;
if SV == 0
sigma2 = out.sigma2;
end

ypredtt = []; % local storage
crps_temp = []; % local storage for crps values

if almonrest == 1
% Get out of sample Almon data
if v < 7
[Xv,~] = midas_dat_r2_final(Xm(1:tin-3+tperiod,xind),grp_idx,poly);
else
    [Xv,~] = midas_dat_r2_final(Xm(1:tin+tperiod,xind),grp_idx,poly);
end
Xv = (Xv);
Xv = Xv(end,:);

else
    if v < 7
    Xv = Xm(tin-3+tperiod,xind);
    else
        Xv = Xm(tin+tperiod,xind);
    end
end

t_cont = 1;
sv_cont = 1;
trend_cont = 0;

% Monte Carlo Integration for predictive distribution
for j = 1:(MCMC)
  
    % Sample g_{t+1}
    if v<7
      g_temp = gout(end,j) + randn*thetaout(j,2);
      tau_temp = tauout(end,j) + exp(0.5*g_temp)*randn;
      g_temp = g_temp + randn*thetaout(j,2);
    else
        g_temp = gout(end,j) + randn*thetaout(j,2);
    end
     % Sample h_{t+1}
     if v<7
      h_temp = hout(end,j) + randn*thetaout(j,1);
      h_temp = h_temp + randn*thetaout(j,1);
     else
         h_temp = hout(end,j) + randn*thetaout(j,1);
     end

     % Sample tau_{t+1}
     if v<7
      tau_temp = tau_temp  + exp(0.5*g_temp)*randn;
     else
         tau_temp = tauout(end,j) + exp(0.5*g_temp)*randn;
     end

     % Sample y_{t+1}

     if t == 1
         t_cont = trnd(nuout(j));
     else
         t_cont = randn;
     end

     if SV == 1
         sv_cont = exp(0.5*h_temp); % add a line for the scale of linear regression
     else
         sv_cont = sqrt(sigma2(j));
     end
     if trend == 1
         trend_cont = tau_temp;
     else 
         trend_cont = out.alpha(j);
     end

      y_temp = trend_cont + Xv*betas_final(j,:)' + sv_cont*t_cont;

      ypredtt= [ypredtt y_temp];
      if t ==1 
      crps_temp = [crps_temp; crps_t(yf(1,1),nuout(j),trend_cont + Xv*betas_final(j,:)',sv_cont)];
      else
          crps_temp = [crps_temp; crps_t(yf(1,1),200,trend_cont + Xv*betas_final(j,:)',sv_cont)];
      end
 end

%%%%%%%%%%%%  Save prediction results by nowcast period %%%%%%%%%%%%%%%
crpsv = [crpsv;mean(crps_temp)];
predv = [predv;ypredtt];
 
 end


%% Save prediction results for each time period
 y_pred = predv;
 pincl(:,:,tperiod) = pincl_temp;


rtresid = [];

for u = 1:vint
test = squeeze(y_pred(u,:,1));
vint1 = mean(test)';
res1 = (vint1 - yf(1,1));
rtresid = [rtresid;res1];
end

crps_all(:,tperiod) = crpsv;
rtresid_all(:,tperiod)= rtresid;
y_pred_all(:,:,tperiod) = y_pred';

end
toc


%%%%%%%% Create Storage Structure after model Estimation and Save
output.resid_all = rtresid_all;
output.crps_all = crps_all;
output.y_pred_all = y_pred_all;
output.d_q = dq_nfor;
output.incl = pincl;
output.yf = yf_all;
output.y = y;


if almonrest == 1
    modelname4 = '_almon';
else
    modelname4 = '_umidas';
end

if ortho_choice == 1
    modelname5 = '_ortho';
else
    modelname5 = '';
end

if group_sparse == 1
    modelname6 = '_groupsparse';
else
    modelname6 = '';
end

if trend == 1
    modelname1 = "_trend";
else
    modelname1 = '';
end

if SV == 1
    modelname2 = "_sv";
else
    modelname2 = '';
end

if t == 1
    modelname3 = "_terr";
else
    modelname3 = '';
end


modname = strcat('results','_iteratednowcasts_newgigg_oldsv_','bg_',num2str(hyperpars(gg,2)),modelname1,modelname2,modelname3,modelname4,modelname5,modelname6,".mat");

save(strcat('Output_iterated','/',modname),"output")

end

delete(gcp('nocreate'))
%% Quick Evaluation
disp(modname)
disp("----------------")
rt_rmsfe1_overnowcasts = std(output.resid_all(:,1:end)')'
rt_rmsfe2_overnowcasts = std(output.resid_all(:,1:51)')'
rt_crps_overnowcasts = mean(crps_all,2)


%%%%%%%%%%%%  Display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format bank
disp('Point evaluation: Average RMSFE across evaluation quarters: by nowcast periods (rows)')
rt_rmsfe1_overnowcasts
disp('Density evaluation: Average CRPS across evaluation quarters: by nowcast periods (rows)')
rt_crps_overnowcasts
disp('Average inclusion probabilities across evaluation quarters, by nowcast periods (row) and indicator (col)')
[names_incl_m; num2cell(mean(pincl,3))]

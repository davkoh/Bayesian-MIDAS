%% Estimation script for "Flexible Bayesian MIDAS: time‑variation, group‑shrinkage and sparsity"
%%% David Kohns, Aalto University
%%% Galina Potjagailo, Bank of England
%%% Date 08.11.2023

%%% This code estimates the T-SVt-BMIDAS model with GIGG prior and ex-post
%%% group sparsification
%   - alternative versions of the model without Trend or SV-t can be
%     selected under "modelling choices"
%   - alternative version without group-sparsification can be selected under "..." 
%   - choice between pseudo real time calendar or simply feeding the latest
%     available data under "Mixed Frequency Estimation Choices"
%%%% ------------------ %%%

clear;
clear all
rng(1,'twister');  %set seed

%% Define directory and upload data
cd '/Users/dk/Documents/GitHub/Bayesian-MIDAS'   %%%%% Specify output directory (replace the current string) %%%%%%%

%mkdir 'Output_v2'
outputfolder = char([cd,'\Outputv2']);  
addpath("Data")
addpath("Matlab")
addpath("R&R")
load("Output/matlab_input.mat")

tperiod = nfor;
v = 1;

for gg = 2 % Loop over all hyperparameters for the GIGG prior
tic
%% Loop over time periods
parfor (tperiod = 1:nfor, cores)
Xf = (Xm(1:tin+tperiod,:));
yf = y(tin+tperiod:end,:); 
T=tin-1+tperiod;

    predv = []; % Locals for the scores
    rtscores = [];
    crpsv = [];
    
    pincl_temp = zeros(vint,max(unique(groupall))); % zeros(vint,size(G,2)); % zeros(vint,size(G,2));
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

Xv = (Xm(1:tin-1+tperiod,xind));

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
input.Y = y(1:tin-1+tperiod,:);
input.X = Xv;
input.burnin = BURNIN;
input.samples = MCMC;
input.btrick = 0;
input.a = repmat(hyperpars(gg,1),sum_grp,1);
input.b = repmat(hyperpars(gg,2),sum_grp,1);
input.standardise = 0;
input.trend = trend;
input.sv = SV;
input.t = t ;

[out] = bmidas_1xdynhs(input);

test = mean(out.tau,2);


% Perform group sparsification
if group_sparse == 1
[beta_out] = group_savs_orth(out.beta,grp_idx_temp');
else
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
hout = out.h;
tauout = out.tau;
nuout = out.nu;

ypredtt = []; % local storage
scores = []; % local storage

if almonrest == 1
% Get out of sample Almon data
[Xv,~] = midas_dat_r2_final(Xm(1:tin+tperiod,xind),grp_idx,poly);
Xv = (Xv);
Xv = Xv(end,:);
else
    Xv = Xm(tin+tperiod,xind);
end

t_cont = 1;
sv_cont = 1;
trend_cont = 1;

% Monte Carlo Integration for predictive distribution
for j = 1:(MCMC)
  
    
     % Sample h_{t+1}
      h_temp = hout(end,j);
     % Sample tau_{t+1}
      tau_temp = tauout(end,j);
     % Sample y_{t+1}

     if t == 1
       t_cont = trnd(nuout(j));
     else 
         t_cont = randn;
     end
     if SV == 1
         sv_cont = exp(0.5*h_temp);
     end
     if trend == 1
         trend_cont = tau_temp;
     end

      y_temp = trend_cont + Xv*betas_final(j,:)' + sv_cont*t_cont;

      ypredtt= [ypredtt y_temp];
  
     
      % Compute Log score here
      test = normpdf(yf,trend_cont+ Xv*betas_final(j,:)',sv_cont*t_cont)/(MCMC);
      scores = [scores;test];

 end

%%%%%%%%%%%%  Save prediction results by nowcast period %%%%%%%%%%%%%%%
crps = pscrps(ypredtt, yf(1,1));
crpsv = [crpsv;crps];
predv = [predv;ypredtt];
rtscores = [rtscores;log(sum(scores))];
 
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
rtlogscores_all(:,tperiod)= rtscores;
y_pred_all(:,:,tperiod) = y_pred';

end
toc


%%%%%%%% Create Storage Structure after model Estimation and Save
output.resid_all = rtresid_all;
output.logscore = rtlogscores_all;
output.crps_all = crps_all;
output.y_pred_all = y_pred_all;
output.d_q = dq_nfor;
output.incl = pincl;
output.yf = yf_all;
output.y = y;


if almonrest == 1
    modelname3 = '_almon';
else
    modelname3 = '_umidas';
end

if ortho_choice == 1
    modelname4 = '_ortho';
else
    modelname4 = '';
end

if group_sparse == 1
    modelname5 = '_groupsparse';
else
    modelname5 = '';
end

modname = strcat('output_hs_sw_corrected','_',num2str(hyperpars(gg,1)),'_',num2str(hyperpars(gg,2)),modelname3,modelname4,modelname5,".mat");

save(strcat(outputfolder,'\',modname),"output")

end

delete(gcp('nocreate'))
%% Quick Evaluation

rt_rmsfe_overnowcasts = std(rtresid_all(6:end-1,1:end)')'
rt_crps_overnowcasts = mean(crps_all(6:end-1,1:36),2)

%%%%%%%%%%%%  Display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format bank
disp('Point evaluation: Average RMSFE across evaluation quarters: by nowcast periods (rows)')
rt_rmsfe_overnowcasts
disp('Density evaluation: Average CRPS across evaluation quarters: by nowcast periods (rows)')
rt_crps_overnowcasts
disp('Average inclusion probabilities across evaluation quarters, by nowcast periods (row) and indicator (col)')
[names_incl_m; num2cell(mean(pincl,3))]

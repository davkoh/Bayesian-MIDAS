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
cd 'D:\Github\Bayesian-MIDAS'   %%%%% Specify output directory (replace the current string) %%%%%%%

%mkdir 'Output'
%outputfolder = char([cd,'\Output']);  
addpath("Data")
addpath("Matlab")

%%%%%%%%%%%%%   Load data   %%%%%%%%%%%%%%%
[data_quarterly, names_q]= xlsread('UK_data_bmidas.xlsx','QuarterlyData','A4:e500');
[data_monthly, names_m]= xlsread('UK_data_bmidas.xlsx','MonthlyData','A4:aq2000');
%%% important: 
%  - excel sheet should be read with variable names 
%  - first data row for data_monthly are transformation indices (will be used in clean_data.m, line 27)

%% ---------- Set Data Choices ---------- %%

%%%%%%%%%%%  Select Sample Period  %%%%%%%%%%%%%%%%%%
beg_s = '31-Dec-1998';    %%% put in first quarter of estimation as "last day - last month of quarter (MMM) - year (YYYY)" 
end_s = '30-Sep-2021';    %%% put in last quarter of estimation as "last day - last month of quarter (MMM) - year (YYYY)" 

beg_eval_per = '31-Mar-2011';  % specify quarter in which to begin evaluation period "last day - last month of quarter (MMM) - year" 
                               % full evaluation can be shut off below via "eval_full==0"
      
%%%%%%%%%%%  select groups of monthly series to include (individual series in each group see below)
sur = 1;     % survey data
act = 1;     % activity and trade data
lab = 1;     % labour market series
pr  = 0;     % prices 
mon = 0;     % money and interest rates
mort = 1;     % mortgages
fin = 0;     % financial indices
ie  = 0;     % inflation expectations
vis = 1;     % VISA consumer spending
                   
%%%%%%%%%%% define data transformations
dyoy     = 0;    % 1- y-o-y growth rates/changes , 0- m-o-m or q-o-q changes
stand    = 0;    % 1- standardise all series, 0- standardise only series in levels

%%%%%%%%%%%  Select Variables  %%%%%%%%%%%%%%%%%%%%%%
%%% define selected monthly series - use name as in first row of names_m!
%   (exclude individual series by copying them out of the brakets and spearating by through three dots)
Var  = {};
Varq = {'GDP_Q'};  ...;'CONS';'HOURS';'INV';
Vars = {'CBI_ES';'CBI_S';'CBI_EO';'PMI_M';'PMI_S';'PMI_C';'GfK'};        % surveys: CBIs,PMIs, GfK
Vara = {'IoP';'IoS';'Exp';'Imp'};                                        % IoP,IoS,Exports,Imports
Varl = {'UR';'EMP';'Vacancies';'Hours'}; ...'AWE';'Claimant'             % UE,EMP,Hours,Vacancies, AWE, Claimant count, 
Varp = {'CPI';'CoreCPI';'RPI';'RPIX';'PPIout';'PPIin';'HPr';'Oil'};      % Prices: CPI,CPI core,RPI,RPIX,PPIout,PPIin,HP,Oil
Varm = {'Money';'BaseRate';'LIBOR';'EXR'};                               % Money: M4,Base rate,LIBOR,Exrate
Varmt = {'Mortgage'};                                                    % Mortgages
Varf = {'FTSE';'FTSE250';'FTSEUK';'SP500';'EuroStoxx';'VIX';'UKVIX'};    % Financial: FTSE all/250/UK,SP500,Euro stoxx, VIX, VIXUK
Vari = {'InflExp5y';'City1y';'City5y'};                                  % Infl expect: 5yr market-based, Citi 1y, City5-10y
Varv = {'VISA'};                                                         % VISA consumer spending

%% ---------- Nowcast evaluation choices ------------------- %%
% Nowcast calendar choice
pseudo_cal = 1;      % 1 = pseudo data release calendar (baseline in paper)
                     % 0 = estimation based on latest available data at time of estimation, 
                     
% Choice of out-of sample evaluation               
eval_full = 1;       % 1 - full evaluation over each quarter in nowcast evaluation period starting in beg_eval_per
                     % 0 - only evaluate over LATEST quarter in the sample      

%% ---------- Stylised Calendar Choices  - if pseudo data release calendar chosen above  ---------- %%
%%% Define publication delays within quarter for the stylised calendar (same structure as for defining the variable names)
% Each variable (same order as defined above) is assigned a publication delay according to the latest month it is available for at the end of a quarter.
% I.e., if variable is available up until June at the end of Q2 it receives a 0, if up until April receives a -2 (June = 0, May=-1, April =-2, March=0, Feb=-1,Jan=-2). 
% Delays should not be lower than amount of months analysed for a nowcast cycles (otherwise no data available).
Var_delay  = [];
Varq_delay = [-2];  
Vars_delay = [0;0;0;-1;-1;-1;0];                                          % surveys: CBIs,PMIs, GfK
Vara_delay = [-2;-2;-2;-2];                                               % IoP,IoS,Exports,Imports
Varl_delay = [-2;-2;-2;-2];                                               % UE,EMP,Hours,Vacancies, AWE, Claimant count,
Varp_delay = [];                                                          % Prices: CPI,CPI core,RPI,RPIX,PPIout,PPIin,HP,Oil                                            
Varm_delay = [];                                                          % Money: M4,Base rate,LIBOR,Exrate
Varmt_delay = [-1];                                                       % Mortgages
Varf_delay = [];                                                          % Financial: FTSE all/250/UK,SP500,Euro stoxx, VIX, VIXUK
Vari_delay = [-1;-1;-1];                                                  % Infl expect: 5yr market-based, Citi 1y, City5-10y
Varv_delay = [-1];                                                        % VISA consumer spending

%%% Define Publication release order for stylised calendar month
% numbering identifies order of publication in an idealised month, same number specified if variables are released on the same release day
Var_pubgroup = [];
Varq_pubgroup = [2];
Vars_pubgroup = [6;6;6;1;1;1;6];                                         % surveys: CBIs,PMIs, GfK
Vara_pubgroup = [3;3;3;3];                                               % IoP,IoS,Exports,Imports
Varl_pubgroup = [4;4;4;4];                                               % UE,EMP,Hours,Vacancies, AWE, Claimant count, 
Varp_pubgroup = [];                                                      % Prices: CPI,CPI core,RPI,RPIX,PPIout,PPIin,HP,Oil
Varm_pubgroup = [];                                                      % Money: M4,Base rate,LIBOR,Exrate
Varmt_pubgroup = [5];                                                    % Mortgages
Varf_pubgroup = [];                                                      % Financial: FTSE all/250/UK,SP500,Euro stoxx, VIX, VIXUK
Vari_pubgroup = [];                                                      % Infl expect: 5yr market-based, Citi 1y, City5-10y
Varv_pubgroup = [5];                                                     % VISA consumer spending          

clearvars input % Ignore this

input.mstart = -3; % Starting month for each nowcast cycle. E.g: choose -3 for start in March if the latest reference month of the quarter is June.
input.mend = 2; % Ending month for each nowcast cycle. E.g: choose 2 for ending nowcasting in August if the reference quarter is June.


%% ---------- Mixed Frequency Estimation Choices ---------- %%
% mixed-frequency lag structure
mismatch = 3; % frequency mismatch between LHS and RHS (3 for quarterly vs monthly)
monthvars = 6; % amount of months related to LHS quarterly variable (multiples of mismatch, max 12)
almonrest = 1; % 1 = use almon lag restrictions (at the moment restricted to a 4th degree with 2 endpoint restrictions), 0 = U-MIDAS
poly = 4; % Polynomial degree for the Almon lag


%% ----------  Estimation and Modeling Choices -------------- %%
% Sampler Info
BURNIN = 5000; % Burnin for Gibbs sampler
MCMC = 5000; % Number of Monte Carlo iterations saved
endpoint = monthvars; % How many monthly lags are related to the LHS

% BMIDAS Choices
trend = 1; % 1 = include trend, 0 = don't include trend
SV = 1; % 1 = include SV, 0 = don't include SV
t = 1; % 1 = include t-errors, 0 = normal errors

% Group-sparsification step                    
group_sparse = 1; % 1 = apply group-sparsity, 
            % 0 = don't apply group-sparsification step

% Parallelisation
cores = 10; % Number of threads for parallelising the nowcast loops (if too high, will default max workers specified by matlab copy)



%%
%%%%%%%%%%% ============================================== %%%%%%%%%%%            
%%%%%%%%%%% ============ AUTOMATIC PART ================== %%%%%%%%%%%
%%%%%%%%%%% ============================================== %%%%%%%%%%%     

%%%%%%%% Does the data cleaning 
clean_data    %%%% transform and plot data, and prepare data for estimation  

%%%%%%% Building pseudo publication Calendar 
%  "input" structure which contains the relevant information in order to construct a pseudo real-time calendar as in the
% paper. Please note, that in the calendar generation code below, it is assumed that the final data publication refers to the quarterly variable coming out.

input.K = size(y_m,2); % number of higher frequency indicators
input.mismatch = mismatch; % mismatch in sampling frequency
input.mlags = monthvars; % number of months used for nowcasting
input.pubdelay = Var_delay; % vector of publication delays of dimension equal to number of higher frequency indicators.
input.pubseq = Var_pubgroup; % vector of groupings that define which variables come out in which order.
[puball groupall] = calendar_gen(input);
                       
%%%% Nowcast calendar definitions: End date Date, Start date and number of forecast periods choice
dqend = d_q(end-1);
dqstart = d_q(1);
mstart = eomdate( dqstart - calmonths(monthvars-1)); % Adjust Monthly series for lags
mend = dqend;
if eval_full ==1
     nfor = size(d_q,1)-find(d_q==beg_eval_per); % Number of nowcast quarters: can be altered by changing the date
else
    nfor =1 ; %% only evaluation for latest quarter
end
tin = size(d_q,1)-1-nfor; % Initial in-sample period 

%% Data Helper (Brings y_m into MIDAS) 
vint = size(puball,1); % number of nowcast periods
missingvalues_mixedfrequency_2 % adjusts the data to starting dates and U-MIDAS sampling

if pseudo_cal ==0
    vint = 1;
    pub_m = avail_ind;
    puball = avail_ind;
end

%%%%%%%%%%%%%%% OUTPUT  Storage

G = size(unique(groupall),2);

% Output Matrices
crps_all = zeros(vint,nfor); % Storage for CRPS values
y_pred_all = zeros(MCMC,vint,nfor); % stores predictive distributions for each nowcast
rtresid_all = zeros(vint,nfor); % stores residuals for each nowcast, based on mean predictive
rtlogscores_all = zeros(vint,nfor); % stores log-scores for each nowcast
yf_all =zeros(nfor,1); % Saves the out-of-sample LHS variable for each quarter
dq_nfor = []; % saves dates of quarters to be nowcasted
pincl = zeros(vint,max(unique(groupall)),nfor);
modall = zeros(vint,MCMC,nfor);

ortho_choice = 1; % group-orthononormalisation, needed for validity of the group-sparsification solution, don't change!
hyperpars = [1/tin,1/tin;1/tin,0.5;1/tin,1;1,1/tin;0.5,1/tin;1,1;0.5,0.5]; % Hyperpriors for the GIGG prior

for gg = 1:size(hyperpars,1) % Loop over all hyperparameters for the GIGG prior
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

[out] = bmidas_svhier(input);


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
gout = out.g;
hout = out.h;
tauout = out.tau;
thetaout = out.theta;
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
trend_cont = 0;

% Monte Carlo Integration for predictive distribution
for j = 1:(MCMC)
  
    % Sample g_{t+1}
      g_temp = gout(end,j) + randn*thetaout(j,2);
     % Sample h_{t+1}
      h_temp = hout(end,j) + randn*thetaout(j,1);
     % Sample tau_{t+1}
      tau_temp = tauout(end,j) + exp(0.5*g_temp)*randn;
     % Sample y_{t+1}

     if t == 1
         t_cont = sqrt((1./gamrnd(nuout(j,1)/2,2/nuout(j,1),1,1)));
     end
     if SV == 1
         sv_cont = exp(0.5*h_temp);
     end
     if trend == 1
         trend_cont = tau_temp;
     end

      y_temp = trend_cont + Xv*betas_final(j,:)' + sv_cont*t_cont*randn;

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

modname = strcat('output','_',num2str(hyperpars(gg,1)),'_',num2str(hyperpars(gg,2)),modelname3,modelname4,modelname5,".mat");

save(strcat(outputfolder,'\',modname),"output")

end

delete(gcp('nocreate'))
%% Quick Evaluation

rt_rmsfe_overnowcasts = std(rtresid_all')'
rt_crps_overnowcasts = mean(crps_all,2)

%%%%%%%%%%%%  Display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format bank
disp('Point evaluation: Average RMSFE across evaluation quarters: by nowcast periods (rows)')
rt_rmsfe_overnowcasts
disp('Density evaluation: Average CRPS across evaluation quarters: by nowcast periods (rows)')
rt_crps_overnowcasts
disp('Average inclusion probabilities across evaluation quarters, by nowcast periods (row) and indicator (col)')
[names_incl_m; num2cell(mean(pincl,3))]

%% Estimation script for "Flexible Bayesian MIDAS: time‑variation, group‑shrinkage and sparsity"
%%% David Kohns, Aalto University
%%% Galina Potjagailo, Bank of England
%%% Date 30.04.2025

%%% This code estimates the T-SVt-BMIDAS model with GIGG prior and ex-post
%%% group sparsification


%%%% ------------------ %%%
clear all
rng(1,'twister');  %set seed

% TODO: add back latest data nowcast, figures only produced when over all
% nowcast periods and quarters

% TODO: the three simple figures for the given model here, and a separate script
% for the paper figures (as previously used)


%% Directories

mkdir 'Output'
outputfolder = char([cd,'\Output']);  
addpath("Data/")
addpath("Matlab/")

%% Define Data

% Data file name
dat_file = 'UK_data_bmidas.xlsx'; %%% put in data file name here

% TODO: bring to newest nowcast version!
% Sample Period
beg_s = '31-Dec-1998';    %%% put in first quarter of estimation as "last day - last month of quarter (MMM) - year (YYYY)" 
end_s = '30-Sep-2022';    %%% put in last quarter of estimation as "last day - last month of quarter (MMM) - year (YYYY)" 

beg_eval_per = '31-Mar-2011';  % specify quarter in which to begin evaluation period "last day - last month of quarter (MMM) - year" 
                               % full evaluation can be shut off below via "eval_full==0"


% Select groups of monthly series to include (individual series in each group see below)
sur = 0;     % survey data: NOT AVAILABLE PUBLICLY
act = 1;     % activity and trade data
lab = 1;     % labour market series
mort = 1;     % mortgages

% Select variables in group to include
%   (exclude individual series by copying them out of the brakets and spearating by through three dots)
Var  = {};
Varq = {'GDP_Q'};  ...;'CONS';'HOURS';'INV';
%Vars = {'CBI_ES';'CBI_S';'CBI_EO';'PMI_M';'PMI_S';'PMI_C';'GfK'};        % surveys: CBIs,PMIs, GfK
Vara = {'IoP';'IoS';'Exp';'Imp'};                                        % IoP,IoS,Exports,Imports
Varl = {'UR';'EMP';'Vacancies';'Hours'}; ...'AWE';'Claimant'             % UE,EMP,Hours,Vacancies, AWE, Claimant count, 
Varmt = {'Mortgage'};                                                    % Mortgages


% Define Transformation
dyoy     = 0;    % 1- y-o-y growth rates/changes , 0- m-o-m or q-o-q changes
stand    = 0;    % 1- standardise all series, 0- standardise only series in levels

% mixed-frequency lag structure
mismatch = 3; % frequency mismatch between LHS and RHS (3 for quarterly vs monthly)
monthvars = 6; % amount of months related to LHS quarterly variable (multiples of mismatch, max 12)
almonrest = 1; % 1 = use almon lag restrictions (at the moment restricted to a 4th degree with 2 endpoint restrictions), 0 = U-MIDAS
poly = 4; % Polynomial degree for the Almon lag

%% Define Stylised Calendar
%%% Define publication delays within quarter for the stylised calendar (same structure as for defining the variable names)
% Each variable (same order as defined above) is assigned a publication delay according to the latest month it is available for at the end of a quarter.
% I.e., for Q2 nowcasts, if variable is available up until June, it receives a 0, if up until April receives a -2 (June = 0, May=-1, April =-2); and for Q1 nowcasts: March=0, Feb=-1,Jan=-2.
Var_delay  = [];
Varq_delay = [-2];  
Vara_delay = [-2;-2;-2;-2];                                               % IoP,IoS,Exports,Imports
Varl_delay = [-2;-2;-2;-2];                                               % UE,EMP,Hours,Vacancies, AWE, Claimant count,
Varmt_delay = [-1];                                                       % Mortgages

%%% Define Publication release order for stylised calendar month
% numbering identifies order of publication in an idealised month, same number specified if variables are released on the same release day
Var_pubgroup = [];
Varq_pubgroup = [2];
Vara_pubgroup = [3;3;3;3];                                               % IoP,IoS,Exports,Imports
Varl_pubgroup = [4;4;4;4];                                               % UE,EMP,Hours,Vacancies, AWE, Claimant count, 
Varmt_pubgroup = [5];                                                    % Mortgages


data_prep

%% MCMC Settings
MCMC =10;
BURNIN = 10;


%% ----- Choose Priors ----- %%

% For all priors below, the default specification follows the T-SV-BMIDAS (GIGG) with
%   1) a hierarchal prior for a_g (see section 2.3),
%   2) fixed hyperparameters for the trend-SV process, 
%   3) fixed hyperparameters for the observation-SV process.

% Any prior choice (if the corresponding model component above was selected) over-rides the
% corresponding default. Any hyper-parameters in the extra hierarchies (e.g. the values for c and d for a_g
% ~ G(c,d)) are set according to the paper.

% Prior choices for exluded model components can be ignored

% BMIDAS (HS) models are defined equivilantly to below with the standard horseshoe
% prior applied to the MIDAS coefficients

%% MIDAS Transformation
    % "almon" = Linear 3rd degree almon polynomial with end-point
    % restrictions
    % "umidas" = linear umidas

midas_transformation_type = "almon";

%% MIDAS Priors

% Choose here which prior is applied to the MIDAS component: "gigg",
% "horseshoe"
    % "gigg" = θ_{k,j} ~ N(0, 𝜗^2𝛾_{k}^2φ_{k,j}^2), γ_{k}^2 ~ G(a_k,1),
    % 𝜑_{k,j}^2 ~ G(b_k,1),
    % "horseshoe" = θ_{k,j} ~ N(0, 𝜗^2φ_{k,j}^2), φ_{k,j} ~ C_+(0,1),
    % "MAL" = following Mogliani & Simoni (2021)

midas_prior = "gigg";


%% GIGG(a_g,b_g) Prior, (a_g,b_g) ~ π()

% Choose hierarchy on a_g and b_g with gigg_type below: "fixed", "hier_ag",
% "hier_bg", "hier_ag_bg"
    % "fixed" = GIGG(a_g = 1/T , b_g = 1/2)
    % "hier_ag" = GIGG(a_g , b_g = 1/2), a_g ~ G(c,d) (ensures relatively
    % high within-group correlation, inference is on group-wise shrinkage)
    % "hier_bg" = GIGG(a_g = 1/T , b_g), b_g ~ G(c,d) (ensures relatively
    % high group-wise shrinkage, inference is on within-group correlation)
    % "hier_ag_bg" = GIGG(a_g , b_g), (a_g,b_g) ~ G(c,d)

% Fixed hyper-parameters (corresponding parameter be over-written if hierichical GIGG is used)
a_g = 1/length(y); % controls goup-level shrinkage
b_g = 0.5; % controls within-group correlation in shrinkage

gigg_type = "fixed"; % Any non-fixed hierarchy selected receives a G(1,2) prior


%% Trend prior, (τ|τ_{t-1}) ~ N(τ_{t-1},σ^{2,τ}_t), σ^{2,τ}_t = exp(g_t), (g_t|g_{t-1},V^2_g) ~ N(g_{t-1},ω^2_g) 
    % Choose which hierarchy on latent trend with trend_type below:
    % "fixed_SV", "PC" , "none"
    % "fixed_SV": ω_g ~ N(0,V^2_{ω_g})  g_0 ~ N(0,V^2_{g_0}), τ_0 ~ N(0,V^2_{τ_0})
    % "none" =  τ_t = α, α ∝ 1
    % "PC": π(V_i^2|ξ_i) for i ∈ {ω_g,g_0,τ_0}

    % For fixed_SV prior
V_omegag = .001;
V_g0 = 0.10;
V_tau0 = 10;

    % For PC prior
xi_g = 0.04; % controls tightness of prior 

trend_type = "fixed_SV";


%% Error process priors, ε^y_t ~ N(0,σ^{2,y}_t) 
    % if sv_obs_type == "none", σ^{2,y}_t = σ^2, σ^2 ∝ 1/σ^2
    % if sv_obs_type == "fixed_SV", σ^{2,y}_t = exp(h_t), (h_t|h_{t-1},V^2_h) ~ N(h_{t-1},ω^2_h) 
    % if sv_obs_type == "PC", σ^{2,y}_t = exp(h_t), (h_t|h_{t-1},V^2_h) ~ N(h_{t-1},ω^2_h)  ,π(V_i^2|ξ_i) for i ∈ {ω_h,h_0}
    % if sv_obs_type == "DHS", σ^{2,y}_t follows a the dynamic horseshoe prior of Kowal et al. (2019)
    % if sv_obs_type == "terr", ε_t ~ t_{ν^y}(0,σ^{2,y}_t) (not available when DHS == 1),
    % prior for ν^y follows recommendations of the paper

% SV Prior, (h|h_{t-1}) ~ N(h_{t-1},σ^{2,h}_t))
    % if "fixed": ω_h ~ N(0,V^2_{ω_h})  h_0 ~ N(0,V^2_{h_0})
V_omegah = .001;
V_h0 = 0.1;

    % if "PC": π(V_i^2|ξ_i) for i ∈ {ω_h,h_0}
xi_h = 0.04; % controls tightness of prior

sv_obs_type = "fixed_SV";

%%%%%%%%%%%%% Automatic from here:

% Output Matrices
crps_all = zeros(vint,nfor); % Storage for CRPS values
y_pred_all = zeros(MCMC,vint,nfor); % stores predictive distributions for each nowcast
rtresid_all = zeros(vint,nfor); % stores residuals for each nowcast, based on mean predictive
rtlogscores_all = zeros(vint,nfor); % stores log-scores for each nowcast
yf_all =zeros(nfor,1); % Saves the out-of-sample LHS variable for each quarter
dq_nfor = []; % saves dates of quarters to be nowcasted
pincl = zeros(vint,max(unique(groupall)),nfor);
modall = zeros(vint,MCMC,nfor);

% Collect prior choices
prior.midas_prior = midas_prior;
prior.gigg_type = gigg_type;
prior.sv_obs_type = sv_obs_type;
prior.a_g = a_g;
prior.b_g = b_g;
prior.trend_type = trend_type;
prior.V_omegag = V_omegag;
prior.V_g0 = V_g0;
prior.V_tau0 =  V_tau0;
prior.xi_g = xi_g;
prior.V_omegah = V_omegah;
prior.V_h0 = V_h0;
prior.xi_h = xi_h;


tic
%% Loop over time periods
for tperiod = 1:nfor

% Update Data and storage locals    
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

% Script that for input into bmidas.m fuction
prep_nowcast_data


%%%%%%%%%%%%%%%%%%  Estimate Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Data
input = [];
input.grp_idx = grp_idx_temp';
if v < 7 % for multiple step ahead forecasts
    input.Y = y(3:tin-2+tperiod,:);
else
    input.Y = y(1:tin-1+tperiod,:);
end
input.X = Xv;
input.burnin = BURNIN;
input.samples = MCMC;
input.btrick = 0;
input.standardise = 0;
input.trend_ind = trend;
input.sv_ind = SV;
input.t_ind = t ;

% Update Prior
input.prior = prior;
input.prior.a_g = repmat(hyperpars(gg,1),sum_grp,1);
input.prior.b_g = repmat(hyperpars(gg,2),sum_grp,1);




[out] = bmidas_wrapper(input,midas_prior);


%%%%%%%%%%%%%%%%%%  Post Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sparsify and find inclusion probabilities
posterior_postprocess

%%%%%%%%%%%%%%%%%%%  Nowcasting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform nowcasting step
nowcast
 
 end


%% Save predictions
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

% TODO: saving option whether they want to save

saveModelOutput(output, ...
                midas_prior, gigg_type, sv_obs_type, ...
                trend_type, midas_transformation_type);







delete(gcp('nocreate'))
%% Quick Evaluation
disp("----------------")
rt_rmsfe1_overnowcasts = std(output.resid_all(:,1:end)')'
rt_rmsfe2_overnowcasts = std(output.resid_all(:,1:51)')'
rt_crps_overnowcasts = mean(crps_all,2)

% TODO: plotting functions here

%%%%%%%%%%%%  Display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format bank
disp('Point evaluation: Average RMSFE across evaluation quarters: by nowcast periods (rows)')
rt_rmsfe1_overnowcasts
disp('Density evaluation: Average CRPS across evaluation quarters: by nowcast periods (rows)')
rt_crps_overnowcasts
disp('Average inclusion probabilities across evaluation quarters, by nowcast periods (row) and indicator (col)')
[names_incl_m; num2cell(mean(pincl,3))]

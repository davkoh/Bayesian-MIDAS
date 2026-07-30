%% Estimation script for "Flexible Bayesian MIDAS: time‑variation, group‑shrinkage and sparsity"
%%% David Kohns, Aalto University
%%% Galina Potjagailo, Bank of England
%%% Date 29.05.2025

%%% This code estimates the T-SVt-BMIDAS model with GIGG prior and ex-post
%%% group sparsification

%%%% ------------------ %%%
clear all
rng(1,'twister');  %set seed




%% Directories

%mkdir 'output'
outputfolder = fullfile(cd,'output');  %"Output"
addpath("Data/")       %"Data/"
addpath("functions/")  %"Matlab/"
save_output=1;

%% Load Data

% Data file name
data_quarterly= readtable('BMIDAS_data.xlsx','Sheet','quarterly','Range','a1:e107');     %BMIDAS_data_restricted.xls - shared publically
data_monthly= readtable('BMIDAS_data.xlsx','Sheet','monthly','Range','a1:s320');          %BMIDAS_data_restricted.xls - shared publically
data_transf = readtable('BMIDAS_data.xlsx','Sheet','monthly_metadata','Range','b1:s2');    %BMIDAS_data_restricted.xls - shared publically


%% DEFINE Sample Period
%%%%% quarters denominated by Mar (Q1), Jun (q2), Sep(q3), Dec (q4)
beg_s = 'Mar-1997';    %%% put in first quarter of estimation as "last day - last month of quarter (MMM) - year (YYYY)" 
end_s = 'Jun-2023';    %%% put in last quarter of estimation as "last day - last month of quarter (MMM) - year (YYYY)" 

beg_eval_per = '31-Mar-2007';  % specify quarter in which to begin evaluation period "last day - last month of quarter (MMM) - year" 
                               % full evaluation can be shut off below via "eval_full==0"

%% Select indicators 

% SELECT groups of monthly series to include (individual series in each group see below)
sur = 1;     % survey data: NOT AVAILABLE PUBLICLY - put to zero if you're using data set without survey data
act = 1;     % activity and trade data
lab = 1;     % labour market series
mort = 1;     % mortgages

% SELECT variables in group to include
%   (exclude individual series by copying them out of the brakets and spearating by through three dots)
Varq = {'GDP_Q'};  ...;'CONS';'HOURS';'INV';
Vars = {'CBI_ES';'CBI_S';'CBI_EO';'PMI_M';'PMI_S';'PMI_C';'GfK'};        % surveys: CBIs,PMIs, GfK
Vara = {'IoP';'IoS';'Exp';'Imp'};                                        % IoP,IoS,Exports,Imports
Varl = {'UR';'EMP';'Vacancies';'Hours'}; ...'AWE';'Claimant'             % UE,EMP,Hours,Vacancies, AWE, Claimant count, 
Varmt = {'Mortgage' };       %                                            % Mortgages

% DEFINE Transformation
dyoy     = 0;    % 1- y-o-y growth rates/changes , 0- m-o-m or q-o-q changes
stand    = 0;    % 1- standardise all series, 0- standardise only series in levels


%% Build data set according to sample selection above
Var=[];
if sur ==1,   Var = [Var; Vars];  end
if act ==1,   Var = [Var; Vara]; end
if lab ==1,   Var = [Var; Varl];  end
if mort ==1,   Var = [Var; Varmt];  end
Var = Var(~cellfun('isempty',Var));
Varq = Varq(~cellfun('isempty',Varq));

% dates and names
d_m = datetime(table2array(data_monthly(:,1)),'InputFormat','MMM-yyyy');
d_m.Format = 'MMM-yyyy';
d_q = datetime(table2array(data_quarterly(:,1)),'InputFormat','MMM-yyyy');
d_q.Format = 'MMM-yyyy';
varnames_qm = [Var',Varq'];
tm_beg=find(ismember(cellstr(d_m),beg_s)); tm_end=find(ismember(cellstr(d_m),end_s));
tq_beg=find(ismember(cellstr(d_q),beg_s)); tq_end=find(ismember(cellstr(d_q),end_s));
d_m=d_m(tm_beg:tm_end);
d_q=d_q(tq_beg:tq_end);

%%% select variables
data_m = table2array(data_monthly(tm_beg:tm_end,Var'));
data_q = table2array(data_quarterly(tq_beg:tq_end,Varq'));
transf = table2array(data_transf(1,Var'));        %transforms for monthly data defined in excel
clear tm_beg tm_end tq_beg tq_end data_monthly data_quarterly Var Varq Varl Vara Varmt Vars 

%%% data transformations applied
[y_m, y_q] = transf_data(data_m, data_q,transf, dyoy, stand);
d_m=d_m(2+11*dyoy:end);    % adjust monthly sample to start at the first quarter without nans
[T,input.K] = size(y_m);
%%%% fill missing values with PCA/ EM algorithm, using 6 principal components
y_m = fillmisspca(y_m, 6) ; 




%% ---------- Stylised release calendar choices ------------------- %%

%%%%%%%%% DEFINE Choice of out-of sample evaluation     %%%%%%%%%%%%%%           
eval_full = 1;       % 1 - full evaluation over each quarter in nowcast evaluation period starting in beg_eval_per
                     % 0 - only evaluate over LATEST quarter in the sample    
%%%%%%%%%%%%%%%%%%%

% DEFINE Nowcast calendar choice
pseudo_cal = 1;      % 1 = pseudo data release calendar (baseline in paper)
                     % 0 = estimation based on latest available data at time of estimation 
                     
if pseudo_cal ==1      % if pseudo realtime calendar - make additional choices          
%%% DEFINE publication delays within quarter for the stylised calendar (same structure as for defining the variable names)
% Each variable (same order as defined above) is assigned a publication delay according to the latest month it is available for at the end of a quarter.
% I.e., for Q2 nowcasts, if variable is available up until June, it receives a 0, if up until April receives a -2 (June = 0, May=-1, April =-2); and for Q1 nowcasts: March=0, Feb=-1,Jan=-2.
    cal.q_delay = [-2];                                                        % quarterly: GDP
    cal.s_delay = [0;0;0;-1;-1;-1;0];                                          % surveys: CBIs,PMIs, GfK
    cal.a_delay = [-2;-2;-2;-2];                                               % IoP,IoS,Exports,Imports
    cal.l_delay = [-2;-2;-2;-2];                                               % UE,EMP,Hours,Vacancies, AWE, Claimant count,
    cal.mt_delay = [-1];                                                       % Mortgages
    input.pubdelay = [cal.q_delay; sur*cal.s_delay; act*cal.a_delay; lab*cal.l_delay; mort*cal.mt_delay];    % combine delays in one vector

%%% DEFINE Publication release order for stylised calendar month
% numbering identifies order of publication in an idealised month, same number specified if variables are released on the same release day
    cal.q_pubseq = [2];
    cal.s_pubseq = [5;5;5;1;1;1;5];                                         % surveys: CBIs,PMIs, GfK
    cal.a_pubseq = [3;3;3;3];                                               % IoP,IoS,Exports,Imports
    cal.l_pubseq = [4;4;4;4];                                               % UE,EMP,Hours,Vacancies, AWE, Claimant count, 
    cal.mt_pubseq = [6];                                                    % Mortgages
    input.pubseq= [cal.q_pubseq; sur*cal.s_pubseq; act*cal.a_pubseq; lab*cal.l_pubseq; mort*cal.mt_pubseq];    % combine pub sequences in one vector

% DEFINE number of monthly lags in quarter    
input.mismatch = 3; % frequency mismatch between LHS and RHS (3 for quarterly vs monthly)
input.mlags = 6; % amount of months related to LHS quarterly variable (multiples of mismatch, max 12) - 3 months for current quarter, and up to three quarters as lagged months
input.qlag=input.mlags/input.mismatch;

% DEFINE span of nowcast calendar
    input.mstart = -2; % Starting month for each nowcast cycle. E.g: choose -3 for start in March if the latest reference month of the quarter is June.
    input.mend = 2; % Ending month for each nowcast cycle. E.g: choose 2 for ending nowcasting in August if the reference quarter is June.

%%%% here calendar for all stylised release dates is being buil
    [puball, groupall] = calendar_gen(input);
    puball = puball(1:end-1,:); % because the calendar adds the last gdp outcome once more to the end
    vint = size(puball,1); % number of nowcast periods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif pseudo_cal==0      %%% here calendar parameters computed when evaluation only for latest nowcast period
    avail_ind = ~isnan(y_m(T-input.mlags+1:T,:));     % Get availability of latest data for nowcasting conditional on availability
    vint = 1;  % number of nowcast periods - here only one (last period) 
    puball = vec(avail_ind)';  %pub_m = avail_ind;       
end

% adjust quarterly data for publication lag
y=y_q(1+input.qlag+3*dyoy:end);  
d_q=d_q(1+input.qlag+3*dyoy:end);   
Tq = size(y);

% Number of nowcast quarters: either only latest or all evaluation period
nfor = 1+ eval_full*(size(d_q,1)-find(d_q==beg_eval_per)-1);
tin = size(d_q,1)-1-nfor; % Initial in-sample period

%%%%%%% monthly data vector brought into UMIDAS form 
Xm = transf_umidas(y,y_m,input);



%%% --------------------------------------------------------------- %%
%%  ---------- Choose MIDAS equation parameters ------------------- %%
%%% --------------------------------------------------------------- %%
midas_type = "almon";   % "almon" = Linear 3rd degree almon polynomial with end-point restrictions
                        % "umidas" = linear umidas
poly = 4; % Polynomial degree for the Almon lag


% Ex-post sparsification choice for the BMIDAS component
post_spars = "yes";
    % yes:                         %%%%%%%%%%%%%%%%%%%% ADD MORE DOCUMENTATION
    % no : ....


%%% ------------------------ %%%
%%  ----- Choose Priors ----- %%
%%% ------------------------- %%

% For all priors below, the default specification follows the T-SV-BMIDAS (GIGG) with
%   1) a hierarchal prior for a_g (see section 2.3),
%   2) fixed hyperparameters for the trend-SV process, 
%   3) fixed hyperparameters for the observation-SV process.

% Any prior choice (if the corresponding model component above was selected) over-rides the
% corresponding default. Any hyper-parameters in the extra hierarchies (e.g. the values for c and d for a_g
% ~ G(c,d)) are set according to the paper.

% Prior choices for excluded model components can be ignored

% BMIDAS (HS) models are defined equivalently to below with the standard horseshoe
% prior applied to the MIDAS coefficients

%% a) MIDAS Priors

% Choose here which prior is applied to the MIDAS component: 3 options
prior.midas = "gigg";
    % "gigg"       :  θ_{k,j} ~ N(0, 𝜗^2𝛾_{k}^2φ_{k,j}^2), γ_{k}^2 ~ G(a_k,1), 𝜑_{k,j}^2 ~ G(b_k,1),
    % "horseshoe"  :  θ_{k,j} ~ N(0, 𝜗^2φ_{k,j}^2), φ_{k,j} ~ C_+(0,1),
    % "MAL"        :  following Mogliani & Simoni (2021)


%% b) GIGG(a_g,b_g) Hierarchy Prior, (a_g,b_g) ~ π()

% Choose among 4 options for hierarchy on a_g and b_g with gigg_type below:
prior.gigg_hyper = "hier_a"; % Any non-fixed hierarchy selected receives a G(1,2) prior
    % "fixed" = GIGG(a_g = 1/T , b_g = 1/2)
    % "hier_a" = GIGG(a_g , b_g = 1/2), a_g ~ G(c,d) (Baseline from paper: ensures relatively
                 % high within-group correlation, inference is on group-wise shrinkage)
    % "hier_b" = GIGG(a_g = 1/T , b_g), b_g ~ G(c,d) (ensures relatively
                  % high group-wise shrinkage, inference is on within-group correlation)
    % "hier_a_b" = GIGG(a_g , b_g), (a_g,b_g) ~ G(c,d)

% If Fixed hyper-parameters (will be over-written if one of hierichical GIGG options is used)
a_g = 1/length(y); % controls goup-level shrinkage
b_g = 0.5; % controls within-group correlation in shrinkage

% if prior is other then GIGG, hyperpara empty
if prior.midas~="gigg", prior.gigg_hyper="";end

%% c) Trend volatility prior, (τ|τ_{t-1}) ~ N(τ_{t-1},σ^{2,τ}_t), σ^{2,τ}_t = exp(g_t), (g_t|g_{t-1},V^2_g) ~ N(g_{t-1},ω^2_g) 

% Choose prior on volatility of latent trend 
prior.trend_sv = "fixed_SV";
    % "none"     :  τ_t = α, α ∝ 1
    % "fixed_SV" :  ω_g ~ N(0,V^2_{ω_g})  g_0 ~ N(0,V^2_{g_0}), τ_0 ~ N(0,V^2_{τ_0})
    % "PC":      :  π(V_i^2|ξ_i) for i ∈ {ω_g,g_0,τ_0}  - Penalised complexity hierarchical prior

    % If fixed_SV prior
V_omegag = .001; 
V_g0 = 10;
V_tau0 = 10;

    % If PC prior
xi_g = 0.04; % controls tightness of prior 

%% d) Stochastic volatility observation equation priors, ε^y_t ~ N(0,σ^{2,y}_t) 

% Choose prior for sotchstic volatility in observation equation
prior.sv_obs = "fixed_SV";
    % "none"     :   σ^{2,y}_t = σ^2, σ^2 ∝ 1/σ^2
    % "fixed_SV" :   σ^{2,y}_t = exp(h_t), (h_t|h_{t-1},V^2_h) ~ N(h_{t-1},ω^2_h) 
    % "PC"       :   σ^{2,y}_t = exp(h_t), (h_t|h_{t-1},V^2_h) ~ N(h_{t-1},ω^2_h)  ,π(V_i^2|ξ_i) for i ∈ {ω_h,h_0} Penalised complexity hierarchical prior
    % "DHS"      :   σ^{2,y}_t follows a the dynamic horseshoe prior of Kowal et al. (2019)

% if "fixed_SV":  (h|h_{t-1}) ~ N(h_{t-1},σ^{2,h}_t));   ω_h ~ N(0,V^2_{ω_h})  h_0 ~ N(0,V^2_{h_0})
V_omegah = .001; 
V_h0 = 10;  

% if "PC": π(V_i^2|ξ_i) for i ∈ {ω_h,h_0}
xi_h = 0.04; % controls tightness of prior


%% e) Tail thickness for errors in observation equation
prior.tail_type = "norm";
    % "norm":  ε_t ~ N_(0,σ^{2,y}_t)
    % "terr":    ε_t ~ t_{ν^y}(0,σ^{2,y}_t) - Student t-errors for fat tails (not available when DHS == 1)


%% MCMC Settings
MCMC =1000;
BURNIN = 1000;


% ----------------------------------------------------------------------------------- % 

% ------------------  Automatic from here  ----------------------------------------- %

% ----------------------------------------------------------------------------------- % 

% Output Matrices
crps_all = zeros(vint,nfor); % Storage for CRPS values
y_pred_all = zeros(MCMC,vint,nfor); % stores predictive distributions for each nowcast
rtresid_all = zeros(vint,nfor); % stores residuals for each nowcast, based on mean predictive
rtlogscores_all = zeros(vint,nfor); % stores log-scores for each nowcast
yf_all =zeros(nfor,1); % Saves the out-of-sample LHS variable for each quarter
dq_nfor = []; % saves dates of quarters to be nowcasted
pincl = zeros(vint,max(unique(groupall)),nfor); % saves the inclusion probabilities
modall = zeros(vint,MCMC,nfor); % saves the model sizes 

% Collect prior coefficients
prior.a_g = a_g;
prior.b_g = b_g;
prior.V_omegag = V_omegag;
prior.V_g0 = V_g0;
prior.V_tau0 =  V_tau0;
prior.xi_g = xi_g;
prior.V_omegah = V_omegah;
prior.V_h0 = V_h0;
prior.xi_h = xi_h;

desiredWorkers = 6;
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', desiredWorkers);
elseif pool.NumWorkers ~= desiredWorkers
    delete(pool);
    parpool('local', desiredWorkers);
end


tic
disp(['Draws for model: ' prior.midas, prior.gigg_hyper,...
    'ex-post sparsification:', post_spars, ', with trend: ' prior.trend_sv,...
    ', observation eq. SV: ', prior.sv_obs, ' error structure: ', prior.tail_type]);

%% Loop over time periods
parfor tperiod = 1:nfor 

% Update Data and storage locals    
Xf = Xm(1:tin+tperiod,:);
yf = y(tin+tperiod:end,:); 
T=tin-1+tperiod;

predv = []; % Locals for the scores
cycpredv = [];
trendv = [];
svv = [];
sv_trendv = [];
rtscores = [];
crpsv = [];
    
pincl_temp = zeros(vint,max(unique(groupall))); 
mod_temp = zeros(vint,MCMC);
            
yf_all(tperiod,:) = yf(1,1);

%% Loop over nowcast periods
for v = 1:vint

    if pseudo_cal == 1
       disp(['Evaluation for time period ' num2str(tperiod) ' of ' num2str(nfor) ', nowcast period ' num2str(v) ' of ' num2str(vint)])    %%% one additional dimension here over which draws are looped?
    else
        disp(['Evaluation for time period ' num2str(tperiod) ' of ' num2str(nfor)])    %%% one additional dimension here over which draws are looped?
    end


% Prepare the nowcast data
nowcast_data = struct('puball',puball,'Xm',Xm,'tperiod',tperiod,'tin',tin,'midas_type',midas_type,'poly',poly,'midas_prior', prior.midas,'v',v,'groupall',groupall);
[Xv,sum_grp,grp_idx_temp,Qj,Lam_inv_sqr,xind,grp_idx] = get_nowcast_data(nowcast_data);


%%%%%%%%%%%%%%%%%%  Estimate Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Data
input = [];
input.grp_idx = grp_idx_temp';
if v < 2 % for multiple step ahead forecasts
    input.Y = y(3:tin-2+tperiod,:);
else
    input.Y = y(1:tin-1+tperiod,:);
end
input.X = Xv;
input.burnin = BURNIN;
input.samples = MCMC;
input.btrick = 0;
input.standardise = 0;
input.prior = prior;
input.prior.a_g = repmat(a_g,sum_grp,1);
input.prior.b_g = repmat(b_g,sum_grp,1);


[out] = bmidas_wrapper(input,prior.midas);


%%%%%%%%%%%%%%%%%%  Post Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sparsify and find inclusion probabilities

post_process_data = struct('post_process',post_spars,'out',out,'grp_idx_temp',grp_idx_temp,...
    'midas_prior',prior.midas,'sum_grp',sum_grp,'Qj',Qj,'Lam_inv_sqr',Lam_inv_sqr,...
    'tin',tin,'v',v,'MCMC',MCMC,'xind',xind,'groupall',groupall,'pincl_temp',pincl_temp,...
    'Xv',Xv);
[betas_final,pincl_temp] = get_sparse_posterior(post_process_data);



%%%%%%%%%%%%%%%%%%%  Nowcasting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform nowcasting step add predv and the other temporaries!! 
nowcast_data = struct('betas_final',betas_final, 'out', out , 'sv_obs_type', prior.sv_obs,...
    'midas_type',midas_type, 'Xm',Xm , 'grp_idx', grp_idx,'poly',poly,'tperiod',tperiod,...
    'xind',xind, 'v',v, 'tail_type',prior.tail_type, 'trend_type',prior.trend_sv,...
    'tin',tin,'MCMC',MCMC,'yf',yf,'crpsv',crpsv,'predv',predv,'midas_prior',prior.midas,...
    'cycpredv',cycpredv, 'trendv',trendv,'svv',svv,'sv_trendv',sv_trendv);
[crpsv, predv,trendv,svv,sv_trendv,cycpredv] = get_nowcasts(nowcast_data);


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

%% Save output
dq_nfor = [dq_nfor d_q(tin+tperiod)];
crps_all(:,tperiod) = crpsv;
rtresid_all(:,tperiod)= rtresid;
y_pred_all(:,:,tperiod) = y_pred';
if ~strcmp(prior.midas,"MAL")
cyc_pred_all(:,:,tperiod) = cycpredv';
if ~strcmp(prior.trend_sv,"none")
    tau_all(:,:,tperiod) = trendv';
    sv_trend_all(:,:,tperiod) = sv_trendv';
end
if ~strcmp(prior.sv_obs,"none")
    sv_all(:,:,tperiod) = svv';
end
end

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
if ~strcmp(prior.trend_sv,"none")
    output.tau_all = tau_all;
    output.sv_trend_all = sv_trend_all;
end
if ~strcmp(prior.sv_obs,"none")
    output.sv_all = sv_all;
end
output.cyc_pred_all = cyc_pred_all;
output.var_names = varnames_qm(1:end-1);

 
% Save Output
output.modelname = strcat('T_',prior.trend_sv, '_OBS_', prior.sv_obs, '_',...
    prior.tail_type, '_BMIDAS_',midas_type,'_', prior.midas,'_', prior.gigg_hyper,...
    '_postspars_',post_spars);

% Create folder name based on model definition
foldername = fullfile(char(outputfolder) , char(output.modelname)) ;
if ~exist(foldername, 'dir'), mkdir(foldername); end % Ensure the folder exists; if not, create it
    
% Save the output structure
if save_output==1
save(fullfile(foldername, strcat('results', '.mat')), 'output');
fprintf('Model output saved in folder: %s, with filename: results.mat\n', foldername);
end


%% Quick Plottiong
close all
%%% Prep: retrieve outputs for desired sub-samples
periods ={'Full', datestr(dq_nfor(1),'mmm-yyyy'), datestr(dq_nfor(end),'mmm-yyyy') ;...
    'GFC', datestr(dq_nfor(1),'mmm-yyyy'),  'Dec-2009' ;...
    'Tranquil', 'Jan-2010', 'Dec-2019';...
    'Pandemic period', 'Jan-2020', datestr(dq_nfor(end),'mmm-yyyy');...
     'Pre pandemic', datestr(dq_nfor(1),'mmm-yyyy'), 'Dec-2019';...
    };
nper = length(periods);
rmsfe = cell(1,nper);
rtcrps = cell(1,nper);
inclprob = cell(1,nper);

% Compute metrics
dq_dates = dateshift(dq_nfor, 'start', 'day');
for i = 1:nper
    period_ind = find(dq_dates >= periods{i,2} & dq_dates <= periods{i,3})';
    rmsfe{i}    = std(output.resid_all(:,period_ind)')';
    rtcrps{i}   = mean(output.crps_all(:,period_ind),2);
    inclprob{i} = mean(output.incl(:,:,period_ind),3);
end

% Figures will be saved into the m


%% Figure 1 - RMSFE and CRPS over time
fig1 = figure;
k=[2;3;4];  % choose which periods to depict - periods(k)

t = tiledlayout(2,length(k), 'TileSpacing','compact', 'Padding','compact');
labels = {'135','120','110','95','85','75','60','50','35','15'};
xticks_vals = [1,3,5,7,9,11,13,15,17,19];

for j=1:length(k) 
    ax1 = nexttile(j); % RMSFE
    plot(ax1, rmsfe{k(j)}, 'k', 'LineWidth', 2);
    if j==1,ylabel(ax1, 'RMSFE', 'FontSize', 16); end
    if j<3,ylim([0;1.5]);end
    set(ax1, 'FontSize', 12)
    title (periods(k(j),1))

    % CRPS
    ax2 = nexttile(j+3);
    plot(ax2, rtcrps{k(j)}, 'k', 'LineWidth', 2);
    if j==1, ylabel(ax2, 'CRPS', 'FontSize', 16); end
    if j<3,ylim([0;1.5]);end
    xlim(ax1, [1;19]); xticks(ax1, xticks_vals); xticklabels(ax1, labels);
    xlim(ax2, [1;19]);  xticks(ax2, xticks_vals); xticklabels(ax2, labels);
    set(ax2, 'FontSize', 12)
end
xlabel(t, 'Days Until GDP Release','FontSize',16);

set(fig1,'PaperOrientation','landscape');
set(fig1,'PaperUnits','normalized');
set(fig1,'PaperPosition', [0 0 1 1]);
% Save
figname = fullfile(foldername, strcat('Figure_Eval_Metrics'));
     print(fig1, figname, '-dpdf');

%% Figure 2 - Heatmap for inclusion probabilities

fig2 = figure;

 YLabels = 1:vint;
    CustomYLabels = string(YLabels);
    CustomYLabels([2 3 5 6 8 9 11 12 14 15 17 18]) = "";
    k=[5;4];  % choose which periods to plot
    for j=1:length(k)
    subplot(2,1,j)
    temp = inclprob{k(j)};
    hm = heatmap(temp,'ColorLimits',[0 1],'CellLabelColor','none');
    colormap(flipud(hot))
    hm.GridVisible = 'off';
    hm.XDisplayLabels = output.var_names;
    hm.YDisplayLabels = CustomYLabels;
    ylabel('Nowcast Periods')
    title(strjoin(periods(k(j),:), ':'));
    hm.FontSize = 15;
    end

set(fig2,'PaperOrientation','landscape');
set(fig2,'PaperUnits','normalized');
set(fig2,'PaperPosition', [0 0 1 1]);

figname = fullfile(foldername, strcat('Figure_InclProb.pdf'));
     print(fig2, figname, '-dpdf');




%% Figure: Trend Decomposition - reads trend, cycle and stochastic volatility from output
% define nowcast periods for which to plot the posterior estimates (between 1 and vint) 
if ~strcmp(prior.trend_sv,"none")
 plot_trend_decomposition(output,[1 12 vint],foldername)
end


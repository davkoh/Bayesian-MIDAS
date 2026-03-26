%% DB Test Tables

%% Load in Model RMSFE and CRPS
clear all

addpath('../Matlab')
% Load Data
load('../Data/UK_dat_2024.mat'); 

%%%%%%%%%%%%%%% Preliminaries + Load models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define number of forecasts and number or nowcast periods
nfor_gfc = 1:12;
nfor_tranq = 13:52;
nfor_covid = 53:65;
nfor_full= 65;
%nfor_pre = 51;  %%% number of nowcast quarter pre-Covid
nper = 19;   %%% number of nowcast periods in each quarter

%%%%%%%%%%%%%%%%  Model names: make sure order corresponds to order in which models are read below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ModNames = ["M&S";"Trend-SV-t-GIGG";"Trend-SV-GIGG";"Trend-GIGG";"SV-t-GIGG";"SV-GIGG";"GIGG";...
    "Trend-SV-GIGG-nospars";"Trend-SV-GIGG-UMIDAS";"Trend-dynHSSV-GIGG";"Trend-dynHSTrend-GIGG";"Trend(PC)-SV(PC)-GIGG";"Trend-SV(PC)-GIGG";"Trend(PC)-SV-GIGG";...
    "Trend-SV-HS";"Trend-SV-SSVS";"Combination";"Trend-SV-Combination";"Trend-SV-PCA"; "Trend(Fixed)-SV-GIGG" ; "HS" ; "SSVS" ;  "PCA-UMIDAS"];
% TODO: String for benchmark model name

nummod =length(ModNames);  % number of models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% TODO: data handling without the unecessary nowcast periods. No change
% could be detected and thrown out, write up once cleanly. Could just
% comment out all but main model, user can do what they want there. 

% TODO: from nq_for, find that subsamples to get rid of hard-coded sub-samp
% definitions.

% TODO: Find index for benchmark model, this should reduce lots of the
% repeated code. 

%%%% load results
mod1 = load('Output_iterated/results_iteratednowcasts_ms2_newdat.mat');        %%%% M&S
mod1= mod1.output; mod1.resid_all([1:5 12 end],:) = []; mod1.crps_all([1:5 12 end],:) = [];
res.mod1=mod1.resid_all; crps.mod1 = mod1.crps_all; wqs.mod1 = calculateWQS(mod1.y_pred_all,mod1.yf,0.05:0.05:0.95,3);wqs.mod1([1:5 12 end],:) =[]; 

mod2 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_trend_sv_terr_almon_ortho_groupsparse.mat') ;    %%% T-SV-t GIGG
mod2= mod2.output; mod2.resid_all([1:5 12 end],:) = []; mod2.crps_all([1:5 12 end],:) = [];
res.mod2=mod2.resid_all; crps.mod2 = mod2.crps_all; wqs.mod2 = calculateWQS(mod2.y_pred_all,mod2.yf,0.05:0.05:0.95,3);wqs.mod2([1:5 12 end],:) =[]; 

mod3 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_trend_sv_almon_ortho_groupsparse.mat');     %%% T-SV GIGG
mod3= mod3.output; mod3.resid_all([1:5 12 end],:) = []; mod3.crps_all([1:5 12 end],:) = [];
res.mod3=mod3.resid_all; crps.mod3 = mod3.crps_all; wqs.mod3 = calculateWQS(mod3.y_pred_all,mod3.yf,0.05:0.05:0.95,3);wqs.mod3([1:5 12 end],:) =[]; 

mod4 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_trend_almon_ortho_groupsparse.mat');     %%% T- GIGG
mod4= mod4.output; mod4.resid_all([1:5 12 end],:) = []; mod4.crps_all([1:5 12 end],:) = [];
res.mod4=mod4.resid_all; crps.mod4 = mod4.crps_all; wqs.mod4 = calculateWQS(mod4.y_pred_all,mod4.yf,0.05:0.05:0.95,3);wqs.mod4([1:5 12 end],:) =[]; 

mod5 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_sv_terr_almon_ortho_groupsparse.mat');     %%% SV-t GIGG
mod5= mod5.output; mod5.resid_all([1:5 12 end],:) = []; mod5.crps_all([1:5 12 end],:) = [];
res.mod5=mod5.resid_all; crps.mod5 = mod5.crps_all; wqs.mod5 = calculateWQS(mod5.y_pred_all,mod5.yf,0.05:0.05:0.95,3);wqs.mod5([1:5 12 end],:) =[]; 

mod6 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_sv_almon_ortho_groupsparse.mat');     %%% SV- GIGG
mod6= mod6.output; mod6.resid_all([1:5 12 end],:) = []; mod6.crps_all([1:5 12 end],:) = [];
res.mod6=mod6.resid_all; crps.mod6 = mod6.crps_all; wqs.mod6 = calculateWQS(mod6.y_pred_all,mod6.yf,0.05:0.05:0.95,3);wqs.mod6([1:5 12 end],:) =[]; 

mod7 = load('Output_iterated/results_iteratednowcasts_oldgigg_plain2_bg_0.5_almon_ortho_groupsparse.mat');     %%% GIGG
mod7= mod7.output; mod7.resid_all([1:5 12 end],:) = []; mod7.crps_all([1:5 12 end],:) = [];
res.mod7=mod7.resid_all; crps.mod7 = mod7.crps_all; wqs.mod7 = calculateWQS(mod7.y_pred_all,mod7.yf,0.05:0.05:0.95,3);wqs.mod7([1:5 12 end],:) =[]; 


mod8 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_trend_sv_almon_ortho.mat');     %%% T-SV GiGG no sparsity
mod8= mod8.output; mod8.resid_all([1:5 12 end],:) = []; mod8.crps_all([1:5 12 end],:) = [];
res.mod8=mod8.resid_all; crps.mod8 = mod8.crps_all; wqs.mod8 = calculateWQS(mod8.y_pred_all,mod8.yf,0.05:0.05:0.95,3);wqs.mod8([1:5 12 end],:) =[]; 

mod9 = load('Output_iterated/results_iteratednowcasts_newgigg_oldsv_bg_0.5_trend_sv_umidas_ortho_groupsparse.mat');     %%% T-SV GiGG UMIDAS 
mod9= mod9.output; mod9.resid_all([1:5 12 end],:) = []; mod9.crps_all([1:5 12 end],:) = [];
res.mod9=mod9.resid_all; crps.mod9 = mod9.crps_all; wqs.mod9 = calculateWQS(mod9.y_pred_all,mod9.yf,0.05:0.05:0.95,3);wqs.mod9([1:5 12 end],:) =[]; 

mod10 = load('Output_iterated/results_iteratednowcasts_gigg_1xdynHSSV_newdat_trend_sv_almon_ortho_groupsparse.mat');     %%% T-GiGG-DynHS in SV
mod10= mod10.output; mod10.resid_all([1:5 12 end],:) = []; mod10.crps_all([1:5 12 end],:) = [];
res.mod10=mod10.resid_all; crps.mod10 = mod10.crps_all; wqs.mod10 = calculateWQS(mod10.y_pred_all,mod10.yf,0.05:0.05:0.95,3);wqs.mod10([1:5 12 end],:) =[]; 

mod11 = load('Output_iterated/results_iteratednowcasts_gigg_1xdynHS_trendnewdat_trend_sv_almon_ortho_groupsparse.mat');     %%% T-GiGG-DynHS in Trend
mod11= mod11.output; mod11.resid_all([1:5 12 end],:) = []; mod11.crps_all([1:5 12 end],:) = [];
res.mod11=mod11.resid_all; crps.mod11 = mod11.crps_all; wqs.mod11 = calculateWQS(mod11.y_pred_all,mod11.yf,0.05:0.05:0.95,3);wqs.mod11([1:5 12 end],:) =[]; 


mod12 = load('Output_iterated/results_gigg_iteratednowcasts_newsv_0.5_trend_sv_almon_ortho_groupsparse.mat');     %%% T-SV-t-GiGG with PC priors in both
mod12= mod12.output; mod12.resid_all([1:5 12 end],:) = []; mod12.crps_all([1:5 12 end],:) = [];
res.mod12=mod12.resid_all; crps.mod12 = mod12.crps_all; wqs.mod12 = calculateWQS(mod12.y_pred_all,mod12.yf,0.05:0.05:0.95,3);wqs.mod12([1:5 12 end],:) =[]; 

mod13 = load('Output_iterated/results_iteratednowcasts_newgigg_newsvSV_bg_0.5_trend_sv_almon_ortho_groupsparse.mat');     %%% T-SV-t-GiGG with PC priors in SV
mod13= mod13.output; mod13.resid_all([1:5 12 end],:) = []; mod13.crps_all([1:5 12 end],:) = [];
res.mod13=mod13.resid_all; crps.mod13 = mod13.crps_all; wqs.mod13 = calculateWQS(mod13.y_pred_all,mod13.yf,0.05:0.05:0.95,3);wqs.mod13([1:5 12 end],:) =[]; 

mod14 = load('Output_iterated/results_iteratednowcasts_gigg_newsvTrend_bg_0.5_trend_sv_almon_ortho_groupsparse.mat');     %%% T-SV-t-GiGG with PC priors in Trend
mod14= mod14.output; mod14.resid_all([1:5 12 end],:) = []; mod14.crps_all([1:5 12 end],:) = [];
res.mod14=mod14.resid_all; crps.mod14 = mod14.crps_all; wqs.mod14 = calculateWQS(mod14.y_pred_all,mod14.yf,0.05:0.05:0.95,3);wqs.mod14([1:5 12 end],:) =[]; 

mod15 = load('Output_iterated/results_iteratednowcasts_hs_oldsv_correcttrend_newdat_trend_sv_almon_groupsparse.mat');     %%% T-SV HS
mod15= mod15.output; mod15.resid_all([1:5 12 end],:) = []; mod15.crps_all([1:5 12 end],:) = [];
res.mod15=mod15.resid_all; crps.mod15 = mod15.crps_all; wqs.mod15 = calculateWQS(mod15.y_pred_all,mod15.yf,0.05:0.05:0.95,3);wqs.mod15([1:5 12 end],:) =[]; 

mod16 = load('Output_iterated/results_iteratednowcasts_ssvs_oldsv_correcttrend_newdat_trend_sv_almon_groupsparse.mat');     %%% T-SV SSVS
mod16= mod16.output; mod16.resid_all([1:5 12 end],:) = []; mod16.crps_all([1:5 12 end],:) = []; 
res.mod16=mod16.resid_all; crps.mod16 = mod16.crps_all; wqs.mod16 = calculateWQS(mod16.y_pred_all,mod16.yf,0.05:0.05:0.95,3); wqs.mod16([1:5 12 end],:) =[]; 

mod17 = load('Output_iterated/results_iteratednowcasts_midascomb_umidas.mat');     %%% Combination plain
mod17= mod17.output; mod17.resid_all([1:5 12 end],:) = []; mod17.crps_all([1:5 12 end],:) = [];
res.mod17=mod17.resid_all; crps.mod17 = mod17.crps_all; wqs.mod17 = calculateWQS(squeeze(mean(mod17.y_pred_all,2)),mod17.yf,0.05:0.05:0.95,3); wqs.mod17([1:5 12 end],:) =[]; 

mod18 = load('Output_iterated/results_iteratednowcasts_midascomb_trend_sv_umidas.mat');     %%% Combination without trend and SV: NEEDS TO BE CHANGED
mod18= mod18.output; mod18.resid_all([1:5 12 end],:) = []; mod18.crps_all([1:5 12 end],:) = [];
res.mod18=mod18.resid_all; crps.mod18 = mod18.crps_all; wqs.mod18 = calculateWQS(squeeze(mean(mod18.y_pred_all,2)),mod18.yf,0.05:0.05:0.95,3); wqs.mod18([1:5 12 end],:) =[]; 

mod19 = load('Output_iterated/results_iteratednowcasts_pca_newsv_newdat_trend_sv_umidas.mat');     %%% PCA Trend-SV-MIDAS
mod19= mod19.output; mod19.resid_all([1:5 12 end],:) = []; mod19.crps_all([1:5 12 end],:) = [];
res.mod19=mod19.resid_all; crps.mod19 = mod19.crps_all; wqs.mod19 = calculateWQS(mod19.y_pred_all,mod19.yf,0.05:0.05:0.95,3); wqs.mod19([1:5 12 end],:) =[]; 

mod20 = load('Output_iterated/results_iteratednowcasts_trendscen_newgigg_oldsv_bg_0.5_trend_sv_almon_ortho_groupsparse.mat');     %%% Trend(fixed)-SV GIGG
mod20= mod20.output; mod20.resid_all([1:5 12 end],:) = []; mod20.crps_all([1:5 12 end],:) = [];
res.mod20=mod20.resid_all; crps.mod20 = mod20.crps_all; wqs.mod20 = calculateWQS(mod20.y_pred_all,mod20.yf,0.05:0.05:0.95,3); wqs.mod20([1:5 12 end],:) =[]; 

mod21 = load('Output_iterated/results_iteratednowcasts_hs_oldsv_correcttrend_newdat_almon_groupsparse.mat');     %%% BMIDAS(HS)
mod21= mod21.output; mod21.resid_all([1:5 12 end],:) = []; mod21.crps_all([1:5 12 end],:) = [];
res.mod21=mod21.resid_all; crps.mod21 = mod21.crps_all; wqs.mod21 = calculateWQS(mod21.y_pred_all,mod21.yf,0.05:0.05:0.95,3); wqs.mod21([1:5 12 end],:) =[]; 

mod22 = load('Output_iterated/results_iteratednowcasts_ssvs_oldsv_correcttrend_newdat_almon_groupsparse.mat');     %%% BMIDAS(SSVS)
mod22= mod22.output; mod22.resid_all([1:5 12 end],:) = []; mod22.crps_all([1:5 12 end],:) = [];
res.mod22=mod22.resid_all; crps.mod22 = mod22.crps_all; wqs.mod22 = calculateWQS(mod22.y_pred_all,mod22.yf,0.05:0.05:0.95,3); wqs.mod22([1:5 12 end],:) =[]; 

mod23 = load('Output_iterated/results_iteratednowcasts_pca_newsv_newdat_umidas.mat');     %%% PCA U-MIDAS
mod23= mod23.output; mod23.resid_all([1:5 12 end],:) = []; mod23.crps_all([1:5 12 end],:) = [];
res.mod23=mod23.resid_all; crps.mod23 = mod23.crps_all; wqs.mod23 = calculateWQS(mod23.y_pred_all,mod23.yf,0.05:0.05:0.95,3); wqs.mod23([1:5 12 end],:) =[]; 




clear mod1 mod2 mod3 mod4 mod5 mod6 mod7 mod8 mod9 mod10 mod11 mod12 mod13 mod14 mod15 mod16 mod17 mod18 mod19 mod20 mod21 mod22 mod23

%% Table 1: Scores Averaged across Nowcast Periods: Against M&S

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================== Retrieve RMSE and CRPS + compute DM =================== %
fm = fieldnames(res);

%RMSE and CRPS, for each nowcast period
for nn=1:numel(fm)
%%%%%%%%%% RMSE and CRPS, average over all nowcast periods
    rmsfe_gfc.(fm{nn}) = std(res.(fm{nn})(:,nfor_gfc),0,2);
    rmsfe_tranq.(fm{nn}) = std(res.(fm{nn})(:,nfor_tranq),0,2);
    rmsfe_covid.(fm{nn}) = std(res.(fm{nn})(:,nfor_covid),0,2);
    rmsfe_full.(fm{nn}) = std(res.(fm{nn})(:,1:nfor_full),0,2);
    rmsfe_pre.(fm{nn}) = std(res.(fm{nn})(:,1:52),0,2);
    crps_gfc.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_gfc),2);
    crps_tranq.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_tranq),2);
    crps_covid.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_covid),2);
    crps_full.(fm{nn}) = mean(crps.(fm{nn})(:,1:nfor_full),2);
    crps_pre.(fm{nn}) = mean(crps.(fm{nn})(:,1:52),2);
    wqs_gfc.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_gfc),2);
    wqs_tranq.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_tranq),2);
    wqs_covid.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_covid),2);
    wqs_full.(fm{nn}) = mean(wqs.(fm{nn})(:,1:nfor_full),2);
    wqs_pre.(fm{nn}) = mean(wqs.(fm{nn})(:,1:52),2);

    for np=1:nper
        %%% over all quarters and for each nowcast period separately
        [dm_pre_each_ms.(fm{nn})(np,:),dmpv_pre_each_ms.(fm{nn})(np,:)]= dmtest_modified(reshape(res.(fm{1})(np,1:52),[],1),reshape(res.(fm{nn})(np,1:52),[],1));
        [dm_pre_each_pca.(fm{nn})(np,:),dmpv_pre_each_pca.(fm{nn})(np,:)]= dmtest_modified(reshape(res.(fm{length(fm)})(np,1:52),[],1),reshape(res.(fm{nn})(np,1:52),[],1));
        [dm_pre_each_ucomb.(fm{nn})(np,:),dmpv_pre_each_ucomb.(fm{nn})(np,:)]= dmtest_modified(reshape(res.(fm{17})(np,1:52),[],1),reshape(res.(fm{nn})(np,1:52),[],1));
        [dm_post_each_ms.(fm{nn})(np,:),dmpv_post_each_ms.(fm{nn})(np,:)]= dmtest_modified(reshape(res.(fm{1})(np,:),[],1),reshape(res.(fm{nn})(np,:),[],1));
        [dm_post_each_pca.(fm{nn})(np,:),dmpv_post_each_pca.(fm{nn})(np,:)]= dmtest_modified(reshape(res.(fm{length(fm)})(np,:),[],1),reshape(res.(fm{nn})(np,:),[],1));
        [dm_post_each_ucomb.(fm{nn})(np,:),dmpv_post_each_ucomb.(fm{nn})(np,:)]= dmtest_modified(reshape(res.(fm{17})(np,:),[],1),reshape(res.(fm{nn})(np,:),[],1));
    end

end

%%%%%%%%%% DM stats Point forecasts
for nn=2:numel(fm)
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_gfc_joint.(fm{nn}),dmpv_gfc_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{1})(:,nfor_gfc),[],1),reshape(res.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_tranq_joint.(fm{nn}),dmpv_tranq_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{1})(:,nfor_tranq),[],1),reshape(res.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_covid_joint.(fm{nn}),dmpv_covid_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{1})(:,nfor_covid),[],1),reshape(res.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_pre_joint.(fm{nn}),dmpv_pre_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{1})(:,1:52),[],1),reshape(res.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_full_joint.(fm{nn}),dmpv_full_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{1})(:,1:nfor_full),[],1),reshape(res.(fm{nn})(:,1:nfor_full),[],1));
end

%%%%%%%%%% DM stats Density forecasts
for nn=2:numel(fm)
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_gfc_joint.(fm{nn}),dmpv_crps_gfc_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{1})(:,nfor_gfc),[],1),reshape(crps.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_tranq_joint.(fm{nn}),dmpv_crps_tranq_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{1})(:,nfor_tranq),[],1),reshape(crps.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_covid_joint.(fm{nn}),dmpv_crps_covid_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{1})(:,nfor_covid),[],1),reshape(crps.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_pre_joint.(fm{nn}),dmpv_crps_pre_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{1})(:,1:52),[],1),reshape(crps.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_full_joint.(fm{nn}),dmpv_crps_full_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{1})(:,1:nfor_full),[],1),reshape(crps.(fm{nn})(:,1:nfor_full),[],1));


    for np=1:nper
        %%% over all quarters pre-pandemic and for each nowcast period separately
        [dm_crps_pre_each_ms.(fm{nn})(np,:),dmpv_crps_pre_each_ms.(fm{nn})(np,:)]= dmtest_modified(reshape(crps.(fm{1})(np,1:52),[],1),reshape(crps.(fm{nn})(np,1:52),[],1));
        [dm_crps_pre_each_pca.(fm{nn})(np,:),dmpv_crps_pre_each_pca.(fm{nn})(np,:)]= dmtest_modified(reshape(crps.(fm{length(fm)})(np,1:52),[],1),reshape(crps.(fm{nn})(np,1:52),[],1));
        [dm_crps_pre_each_ucomb.(fm{nn})(np,:),dmpv_crps_pre_each_ucomb.(fm{nn})(np,:)]= dmtest_modified(reshape(crps.(fm{17})(np,1:52),[],1),reshape(crps.(fm{nn})(np,1:52),[],1));
        %%% over all quarters and for each nowcast period separately
        [dm_crps_post_each_ms.(fm{nn})(np,:),dmpv_crps_post_each_ms.(fm{nn})(np,:)]= dmtest_modified(reshape(crps.(fm{1})(np,:),[],1),reshape(crps.(fm{nn})(np,:),[],1));
        [dm_crps_post_each_pca.(fm{nn})(np,:),dmpv_crps_post_each_pca.(fm{nn})(np,:)]= dmtest_modified(reshape(crps.(fm{length(fm)})(np,:),[],1),reshape(crps.(fm{nn})(np,:),[],1));
        [dm_crps_post_each_ucomb.(fm{nn})(np,:),dmpv_crps_post_each_ucomb.(fm{nn})(np,:)]= dmtest_modified(reshape(crps.(fm{17})(np,:),[],1),reshape(crps.(fm{nn})(np,:),[],1));
    end
    
end

%%%%%%%%%% DM stats Quantile forecasts
for nn=2:numel(fm)
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_gfc_joint.(fm{nn}),dmpv_wqs_gfc_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{1})(:,nfor_gfc),[],1),reshape(wqs.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_tranq_joint.(fm{nn}),dmpv_wqs_tranq_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{1})(:,nfor_tranq),[],1),reshape(wqs.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_covid_joint.(fm{nn}),dmpv_wqs_covid_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{1})(:,nfor_covid),[],1),reshape(wqs.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_pre_joint.(fm{nn}),dmpv_wqs_pre_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{1})(:,1:52),[],1),reshape(wqs.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_full_joint.(fm{nn}),dmpv_wqs_full_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{1})(:,1:nfor_full),[],1),reshape(wqs.(fm{nn})(:,1:nfor_full),[],1));


    for np=1:nper
        %%% over all quarters pre-pandemic and for each nowcast period separately
        [dm_wqs_pre_each_ms.(fm{nn})(np,:),dmpv_wqs_pre_each_ms.(fm{nn})(np,:)]= dmtest_modified(reshape(wqs.(fm{1})(np,1:52),[],1),reshape(wqs.(fm{nn})(np,1:52),[],1));
        [dm_wqs_pre_each_pca.(fm{nn})(np,:),dmpv_wqs_pre_each_pca.(fm{nn})(np,:)]= dmtest_modified(reshape(wqs.(fm{length(fm)})(np,1:52),[],1),reshape(wqs.(fm{nn})(np,1:52),[],1));
        [dm_wqs_pre_each_ucomb.(fm{nn})(np,:),dmpv_wqs_pre_each_ucomb.(fm{nn})(np,:)]= dmtest_modified(reshape(wqs.(fm{17})(np,1:52),[],1),reshape(wqs.(fm{nn})(np,1:52),[],1));
        %%% over all quarters and for each nowcast period separately
        [dm_wqs_post_each_ms.(fm{nn})(np,:),dmpv_wqs_post_each_ms.(fm{nn})(np,:)]= dmtest_modified(reshape(wqs.(fm{1})(np,:),[],1),reshape(wqs.(fm{nn})(np,:),[],1));
        [dm_wqs_post_each_pca.(fm{nn})(np,:),dmpv_wqs_post_each_pca.(fm{nn})(np,:)]= dmtest_modified(reshape(wqs.(fm{length(fm)})(np,:),[],1),reshape(wqs.(fm{nn})(np,:),[],1));
        [dm_wqs_post_each_ucomb.(fm{nn})(np,:),dmpv_wqs_post_each_ucomb.(fm{nn})(np,:)]= dmtest_modified(reshape(wqs.(fm{17})(np,:),[],1),reshape(wqs.(fm{nn})(np,:),[],1));
    end
    
end



temp = [];
for nn=2:numel(fm)   
        temp =[temp; [mean(rmsfe_gfc.(fm{nn}))/mean(rmsfe_gfc.(fm{1})),mean(rmsfe_tranq.(fm{nn}))/mean(rmsfe_tranq.(fm{1})),mean(rmsfe_covid.(fm{nn}))/mean(rmsfe_covid.(fm{1})),mean(rmsfe_pre.(fm{nn}))/mean(rmsfe_pre.(fm{1})), mean(rmsfe_full.(fm{nn}))/mean(rmsfe_full.(fm{1})),...
        mean(crps_gfc.(fm{nn}))/mean(crps_gfc.(fm{1})), mean(crps_tranq.(fm{nn}))/mean(crps_tranq.(fm{1})),mean(crps_covid.(fm{nn}))/mean(crps_covid.(fm{1})),mean(crps_pre.(fm{nn}))/mean(crps_pre.(fm{1})),mean(crps_full.(fm{nn}))/mean(crps_full.(fm{1})),...
        mean(wqs_gfc.(fm{nn}))/mean(wqs_gfc.(fm{1})), mean(wqs_tranq.(fm{nn}))/mean(wqs_tranq.(fm{1})),mean(wqs_covid.(fm{nn}))/mean(wqs_covid.(fm{1})),mean(wqs_pre.(fm{nn}))/mean(wqs_pre.(fm{1})),mean(wqs_full.(fm{nn}))/mean(wqs_full.(fm{1}))],...
        dm_gfc_joint.(fm{nn}), dmpv_gfc_joint.(fm{nn}), ...
        dm_tranq_joint.(fm{nn}), dmpv_tranq_joint.(fm{nn}), ...
        dm_covid_joint.(fm{nn}), dmpv_covid_joint.(fm{nn}), ...
        dm_pre_joint.(fm{nn}), dmpv_pre_joint.(fm{nn}), ...
        dm_full_joint.(fm{nn}), dmpv_full_joint.(fm{nn}), ...
        dm_crps_gfc_joint.(fm{nn}), dmpv_crps_gfc_joint.(fm{nn}),...
        dm_crps_tranq_joint.(fm{nn}), dmpv_crps_tranq_joint.(fm{nn}),...
        dm_crps_covid_joint.(fm{nn}), dmpv_crps_covid_joint.(fm{nn}),...
        dm_crps_pre_joint.(fm{nn}), dmpv_crps_pre_joint.(fm{nn}),...
        dm_crps_full_joint.(fm{nn}), dmpv_crps_full_joint.(fm{nn}),...
        dm_wqs_gfc_joint.(fm{nn}), dmpv_wqs_gfc_joint.(fm{nn}),...
        dm_wqs_tranq_joint.(fm{nn}), dmpv_wqs_tranq_joint.(fm{nn}),...
        dm_wqs_covid_joint.(fm{nn}), dmpv_wqs_covid_joint.(fm{nn}),...
        dm_wqs_pre_joint.(fm{nn}), dmpv_wqs_pre_joint.(fm{nn}),...
        dm_wqs_full_joint.(fm{nn}), dmpv_wqs_full_joint.(fm{nn})];
end


restable = [[ mean(rmsfe_gfc.(fm{1})), mean(rmsfe_tranq.(fm{1})),mean(rmsfe_covid.(fm{1})), mean(rmsfe_pre.(fm{1})) , mean(rmsfe_full.(fm{1})),mean(crps_gfc.(fm{1})), mean(crps_tranq.(fm{1})),mean(crps_covid.(fm{1})), mean(crps_pre.(fm{1})) ,mean(crps_full.(fm{1})),mean(wqs_gfc.(fm{1})), mean(wqs_tranq.(fm{1})),mean(wqs_covid.(fm{1})), mean(wqs_pre.(fm{1})) ,mean(wqs_full.(fm{1})),nan(1,30)]; temp];

Scores1 = ["RMSFE_GFC";"RMSFE_Tranquil";"RMSFE_Covid";"RMSFE_Pre";"RMSFE_Full";"CRPS_GFC";"CRPS_Tranquil";"CRPS_Covid";"CRPS_Pre";"CRPS_Full";"WQS_GFC";"WQS_Tranquil";"WQS_Covid";"WQS_Pre";"WQS_Full";"DM_joint_GFC";"DM_pv_joint_GFC";"DM_joint_Tranquil";"DM_pv_joint_Tranquil";"DM_joint_Covid";"DM_pv_joint_Covid";"DM_joint_Pre";"DM_pv_joint_Pre";"DM_joint_Full";"DM_pv_joint_Full";"DM_CRPS_joint_GFC";"DM_CRPS_pv_joint_GFC";"DM_CRPS_joint_Tranquil";"DM_CRPS_pv_joint_Tranquil";"DM_CRPS_joint_Covid";"DM_CRPS_pv_joint_Covid";"DM_CRPS_joint_Pre";"DM_CRPS_pv_joint_Pre";"DM_CRPS_joint_Full";"DM_CRPS_pv_joint_Full";"DM_WQS_joint_GFC";"DM_WQS_pv_joint_GFC";"DM_WQS_joint_Tranquil";"DM_WQS_pv_joint_Tranquil";"DM_WQS_joint_Covid";"DM_WQS_pv_joint_Covid";"DM_WQS_joint_Pre";"DM_WQS_pv_joint_Pre";"DM_WQS_joint_Full";"DM_WQS_pv_joint_Full"];
tab1 = array2table(restable,'VariableNames',Scores1,'RowNames',ModNames);


%% Table 2: Scores Per Nowcast Period
Nowcast_periods =string((1:nper)');
%%%% Point forecasts
dm_pre_nowcast_ms= [];
dm_pre_nowcast_pca= [];
dm_pre_nowcast_ucomb= [];
dmpv_pre_nowcast_ms = [];
dmpv_pre_nowcast_pca = [];
dmpv_pre_nowcast_ucomb = [];
dm_post_nowcast_ms = [];
dm_post_nowcast_pca = [];
dm_post_nowcast_ucomb = [];
dmpv_post_nowcast_ms = [];
dmpv_post_nowcast_pca = [];
dmpv_post_nowcast_ucomb = [];
rmse_each_pre_ms=nan(nummod-1,nper);
rmse_each_pre_pca=nan(nummod-1,nper);
rmse_each_pre_ucomb=nan(nummod-1,nper);
rmse_each_post_ms=nan(nummod-1,nper);
rmse_each_post_pca=nan(nummod-1,nper);
rmse_each_post_ucomb=nan(nummod-1,nper);
for nn=2:numel(fm)
    dm_pre_nowcast_ms = [dm_pre_nowcast_ms;dm_pre_each_ms.(fm{nn})'];
    dm_pre_nowcast_pca = [dm_pre_nowcast_pca;dm_pre_each_pca.(fm{nn})'];
    dm_pre_nowcast_ucomb = [dm_pre_nowcast_ucomb;dm_pre_each_ucomb.(fm{nn})'];
    dmpv_pre_nowcast_ms = [dmpv_pre_nowcast_ms;dmpv_pre_each_ms.(fm{nn})'];
    dmpv_pre_nowcast_pca = [dmpv_pre_nowcast_pca;dmpv_pre_each_pca.(fm{nn})'];
    dmpv_pre_nowcast_ucomb = [dmpv_pre_nowcast_ucomb;dmpv_pre_each_ucomb.(fm{nn})'];
    dm_post_nowcast_ms = [dm_post_nowcast_ms;dm_post_each_ms.(fm{nn})'];
    dm_post_nowcast_pca = [dm_post_nowcast_pca;dm_post_each_pca.(fm{nn})'];
    dm_post_nowcast_ucomb = [dm_post_nowcast_ucomb;dm_post_each_ucomb.(fm{nn})'];
    dmpv_post_nowcast_ms = [dmpv_post_nowcast_ms;dmpv_post_each_ms.(fm{nn})'];
    dmpv_post_nowcast_pca = [dmpv_post_nowcast_pca;dmpv_post_each_pca.(fm{nn})'];
    dmpv_post_nowcast_ucomb = [dmpv_post_nowcast_ucomb;dmpv_post_each_ucomb.(fm{nn})'];
    for np=1:nper
    rmse_each_pre_ms(nn-1,np) = rmsfe_pre.(fm{nn})(np)/rmsfe_pre.(fm{1})(np);
    rmse_each_pre_pca(nn-1,np) = rmsfe_pre.(fm{nn})(np)/rmsfe_pre.(fm{length(fm)})(np);
    rmse_each_pre_ucomb(nn-1,np) = rmsfe_pre.(fm{nn})(np)/rmsfe_pre.(fm{17})(np);
    rmse_each_full_ms(nn-1,np) = rmsfe_full.(fm{nn})(np)/rmsfe_full.(fm{1})(np);
    rmse_each_full_pca(nn-1,np) = rmsfe_full.(fm{nn})(np)/rmsfe_full.(fm{length(fm)})(np);
    rmse_each_full_ucomb(nn-1,np) = rmsfe_full.(fm{nn})(np)/rmsfe_full.(fm{17})(np);
    end
end
tab_rmse_each_pre_ms = array2table(rmse_each_pre_ms,'RowNames',ModNames(2:end));
tab_rmse_each_pre_pca = array2table(rmse_each_pre_pca,'RowNames',ModNames(2:end));
tab_rmse_each_pre_ucomb = array2table(rmse_each_pre_ucomb,'RowNames',ModNames(2:end));
tab_rmse_each_full_ms = array2table(rmse_each_full_ms,'RowNames',ModNames(2:end));
tab_rmse_each_full_pca = array2table(rmse_each_full_pca,'RowNames',ModNames(2:end));
tab_rmse_each_full_ucomb= array2table(rmse_each_full_ucomb,'RowNames',ModNames(2:end));
tab_dm_each_pre_ms = array2table(dm_pre_nowcast_ms,'RowNames',ModNames(2:end));
tab_dm_each_pre_pca = array2table(dm_pre_nowcast_pca,'RowNames',ModNames(2:end));
tab_dm_each_pre_ucomb= array2table(dm_pre_nowcast_ucomb,'RowNames',ModNames(2:end));
tab_dm_each_full_ms = array2table(dm_post_nowcast_ms,'RowNames',ModNames(2:end));
tab_dm_each_full_pca = array2table(dm_post_nowcast_pca,'RowNames',ModNames(2:end));
tab_dm_each_full_ucomb= array2table(dm_post_nowcast_ucomb,'RowNames',ModNames(2:end));
tab_dmpv_each_pre_ms = array2table(dmpv_pre_nowcast_ms,'RowNames',ModNames(2:end));
tab_dmpv_each_pre_pca = array2table(dmpv_pre_nowcast_pca,'RowNames',ModNames(2:end));
tab_dmpv_each_pre_ucomb= array2table(dmpv_pre_nowcast_ucomb,'RowNames',ModNames(2:end));
tab_dmpv_each_full_ms = array2table(dmpv_post_nowcast_ms,'RowNames',ModNames(2:end));
tab_dmpv_each_full_pca = array2table(dmpv_post_nowcast_pca,'RowNames',ModNames(2:end));
tab_dmpv_each_full_ucomb= array2table(dmpv_post_nowcast_ucomb,'RowNames',ModNames(2:end));

% Tabs for Periods RMSE, CRPS and WQS Score for T-SV-GIGG, PCA and MS for periods
tab_rmsfe_full = array2table([rmsfe_full.mod3 rmsfe_full.mod1 rmsfe_full.mod20 rmsfe_full.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_rmsfe_gfc = array2table([rmsfe_gfc.mod3 rmsfe_gfc.mod1 rmsfe_gfc.mod20 rmsfe_gfc.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_rmsfe_tranq = array2table([rmsfe_tranq.mod3 rmsfe_tranq.mod1 rmsfe_tranq.mod20 rmsfe_tranq.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_rmsfe_covid= array2table([rmsfe_covid.mod3 rmsfe_covid.mod1 rmsfe_covid.mod20 rmsfe_covid.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_rmsfe_pre = array2table([rmsfe_pre.mod3 rmsfe_pre.mod1 rmsfe_pre.mod20 rmsfe_pre.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);

tab_crps_full = array2table([crps_full.mod3 crps_full.mod1 crps_full.mod20 crps_full.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_crps_gfc = array2table([crps_gfc.mod3 crps_gfc.mod1 crps_gfc.mod20 crps_gfc.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_crps_tranq = array2table([crps_tranq.mod3 crps_tranq.mod1 crps_tranq.mod20 crps_tranq.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_crps_covid= array2table([crps_covid.mod3 crps_covid.mod1 crps_covid.mod20 crps_covid.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_crps_pre = array2table([crps_pre.mod3 crps_pre.mod1 crps_pre.mod20 crps_pre.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);

tab_wqs_full = array2table([wqs_full.mod3 wqs_full.mod1 wqs_full.mod20 wqs_full.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_wqs_gfc = array2table([wqs_gfc.mod3 wqs_gfc.mod1 wqs_gfc.mod20 wqs_gfc.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_wqs_tranq = array2table([wqs_tranq.mod3 wqs_tranq.mod1 wqs_tranq.mod20 wqs_tranq.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_wqs_covid= array2table([wqs_covid.mod3 wqs_covid.mod1 wqs_covid.mod20 wqs_covid.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);
tab_wqs_pre = array2table([wqs_pre.mod3 wqs_pre.mod1 wqs_pre.mod20 wqs_pre.mod17],"VariableNames",[ModNames(3);ModNames(1);ModNames(length(fm));ModNames(17)]);



% Tabs for Average by Period for pre Covid RMSE and 
tab_rmse_average = [mean(rmsfe_gfc.(fm{1})) mean(rmsfe_tranq.(fm{1})) mean(rmsfe_covid.(fm{1})) mean(rmsfe_pre.(fm{1})) mean(rmsfe_full.(fm{1})) mean(crps_gfc.(fm{1}))  mean(crps_tranq.(fm{1})) mean(crps_covid.(fm{1})) mean(crps_pre.(fm{1})) mean(crps_full.(fm{1})) mean(wqs_gfc.(fm{1}))  mean(wqs_tranq.(fm{1})) mean(wqs_covid.(fm{1})) mean(wqs_pre.(fm{1})) mean(wqs_full.(fm{1}));...
    mean(rmsfe_gfc.(fm{3}))  mean(rmsfe_tranq.(fm{3})) mean(rmsfe_covid.(fm{3})) mean(rmsfe_pre.(fm{3}))  mean(rmsfe_full.(fm{3})) mean(crps_gfc.(fm{3}))  mean(crps_tranq.(fm{3})) mean(crps_covid.(fm{3})) mean(crps_pre.(fm{3})) mean(crps_full.(fm{3})) mean(wqs_gfc.(fm{3}))  mean(wqs_tranq.(fm{3})) mean(wqs_covid.(fm{3})) mean(wqs_pre.(fm{3})) mean(wqs_full.(fm{3}));...
    mean(rmsfe_gfc.(fm{length(fm)}))  mean(rmsfe_tranq.(fm{length(fm)})) mean(rmsfe_covid.(fm{length(fm)})) mean(rmsfe_pre.(fm{length(fm)}))  mean(rmsfe_full.(fm{length(fm)})) mean(crps_gfc.(fm{length(fm)}))  mean(crps_tranq.(fm{length(fm)})) mean(crps_covid.(fm{length(fm)})) mean(crps_pre.(fm{length(fm)})) mean(crps_full.(fm{length(fm)}))  mean(wqs_gfc.(fm{length(fm)}))  mean(wqs_tranq.(fm{length(fm)})) mean(wqs_covid.(fm{length(fm)})) mean(wqs_pre.(fm{length(fm)})) mean(wqs_full.(fm{length(fm)}));...
    mean(rmsfe_gfc.(fm{17}))  mean(rmsfe_tranq.(fm{17})) mean(rmsfe_covid.(fm{17})) mean(rmsfe_pre.(fm{17}))  mean(rmsfe_full.(fm{17})) mean(crps_gfc.(fm{17}))  mean(crps_tranq.(fm{17})) mean(crps_covid.(fm{17})) mean(crps_pre.(fm{17})) mean(crps_full.(fm{17})) mean(wqs_gfc.(fm{17}))  mean(wqs_tranq.(fm{17})) mean(wqs_covid.(fm{17})) mean(wqs_pre.(fm{17})) mean(wqs_full.(fm{17}))];

Scores1 = ["RMSFE_GFC";"RMSFE_Tranquil";"RMSFE_Covid";"RMSFE_Pre";"RMSFE_Full";"CRPS_GFC";"CRPS_Tranquil";"CRPS_Covid";"CRPS_Pre";"CRPS_Full";"WQS_GFC";"WQS_Tranquil";"WQS_Covid";"WQS_Pre";"WQS_Full"];
tab3 = array2table(tab_rmse_average,'VariableNames',Scores1,'RowNames',[ModNames(1);ModNames(3);ModNames(length(fm));ModNames(17)]);


%%%% Density forecasts
dm_crps_pre_nowcast_ms= [];
dm_crps_pre_nowcast_pca= [];
dm_crps_pre_nowcast_ucomb= [];
dmpv_crps_pre_nowcast_ms = [];
dmpv_crps_pre_nowcast_pca = [];
dmpv_crps_pre_nowcast_ucomb = [];
dm_crps_post_nowcast_ms = [];
dm_crps_post_nowcast_pca = [];
dm_crps_post_nowcast_ucomb = [];
dmpv_crps_post_nowcast_ms = [];
dmpv_crps_post_nowcast_pca = [];
dmpv_crps_post_nowcast_ucomb = [];
crps_each_pre_ms=nan(nummod-1,nper);
crps_each_pre_pca=nan(nummod-1,nper);
crps_each_pre_ucomb=nan(nummod-1,nper);
crps_each_post_ms=nan(nummod-1,nper);
crps_each_post_pca=nan(nummod-1,nper);
crps_each_post_ucomb=nan(nummod-1,nper);
for nn=2:numel(fm)
    dm_crps_pre_nowcast_ms = [dm_crps_pre_nowcast_ms;dm_crps_pre_each_ms.(fm{nn})'];
    dm_crps_pre_nowcast_pca = [dm_crps_pre_nowcast_pca;dm_crps_pre_each_pca.(fm{nn})'];
    dm_crps_pre_nowcast_ucomb = [dm_crps_pre_nowcast_ucomb;dm_crps_pre_each_ucomb.(fm{nn})'];
    dmpv_crps_pre_nowcast_ms = [dmpv_crps_pre_nowcast_ms;dmpv_crps_pre_each_ms.(fm{nn})'];
    dmpv_crps_pre_nowcast_pca = [dmpv_crps_pre_nowcast_pca;dmpv_crps_pre_each_pca.(fm{nn})'];
    dmpv_crps_pre_nowcast_ucomb = [dmpv_crps_pre_nowcast_ucomb;dmpv_crps_pre_each_ucomb.(fm{nn})'];
    dm_crps_post_nowcast_ms = [dm_crps_post_nowcast_ms;dm_crps_post_each_ms.(fm{nn})'];
    dm_crps_post_nowcast_pca = [dm_crps_post_nowcast_pca;dm_crps_post_each_pca.(fm{nn})'];
    dm_crps_post_nowcast_ucomb = [dm_crps_post_nowcast_ucomb;dm_crps_post_each_ucomb.(fm{nn})'];
    dmpv_crps_post_nowcast_ms = [dmpv_crps_post_nowcast_ms;dmpv_crps_post_each_ms.(fm{nn})'];
    dmpv_crps_post_nowcast_pca = [dmpv_crps_post_nowcast_pca;dmpv_crps_post_each_pca.(fm{nn})'];
    dmpv_crps_post_nowcast_ucomb= [dmpv_crps_post_nowcast_ucomb;dmpv_crps_post_each_ucomb.(fm{nn})'];
    for np=1:nper
    crps_each_pre_ms(nn-1,np) = crps_pre.(fm{nn})(np)/crps_pre.(fm{1})(np);
    crps_each_pre_pca(nn-1,np) = crps_pre.(fm{nn})(np)/crps_pre.(fm{length(fm)})(np);
    crps_each_pre_ucomb(nn-1,np) = crps_pre.(fm{nn})(np)/crps_pre.(fm{17})(np);
    crps_each_full_ms(nn-1,np) = crps_full.(fm{nn})(np)/crps_full.(fm{1})(np);
    crps_each_full_pca(nn-1,np) = crps_full.(fm{nn})(np)/crps_full.(fm{length(fm)})(np);
    crps_each_full_ucomb(nn-1,np) = crps_full.(fm{nn})(np)/crps_full.(fm{17})(np);
    end
end
tab_crps_each_pre_ms = array2table(crps_each_pre_ms,'RowNames',ModNames(2:end));
tab_crps_each_pre_pca = array2table(crps_each_pre_pca,'RowNames',ModNames(2:end));
tab_crps_each_pre_ucomb = array2table(crps_each_pre_ucomb,'RowNames',ModNames(2:end));
tab_crps_each_full_ms = array2table(crps_each_full_ms,'RowNames',ModNames(2:end));
tab_crps_each_full_pca = array2table(crps_each_full_pca,'RowNames',ModNames(2:end));
tab_crps_each_full_ucomb = array2table(crps_each_full_ucomb,'RowNames',ModNames(2:end));
tab_crps_dm_each_pre_ms = array2table(dm_crps_pre_nowcast_ms,'RowNames',ModNames(2:end));
tab_crps_dm_each_pre_pca = array2table(dm_crps_pre_nowcast_pca,'RowNames',ModNames(2:end));
tab_crps_dm_each_pre_ucomb= array2table(dm_crps_pre_nowcast_ucomb,'RowNames',ModNames(2:end));
tab_crps_dm_each_full_ms = array2table(dm_crps_post_nowcast_ms,'RowNames',ModNames(2:end));
tab_crps_dm_each_full_pca = array2table(dm_crps_post_nowcast_pca,'RowNames',ModNames(2:end));
tab_crps_dm_each_full_ucomb = array2table(dm_crps_post_nowcast_ucomb,'RowNames',ModNames(2:end));
tab_crps_dmpv_each_pre_ms = array2table(dmpv_crps_pre_nowcast_ms,'RowNames',ModNames(2:end));
tab_crps_dmpv_each_pre_pca = array2table(dmpv_crps_pre_nowcast_pca,'RowNames',ModNames(2:end));
tab_crps_dmpv_each_pre_ucomb = array2table(dmpv_crps_pre_nowcast_ucomb,'RowNames',ModNames(2:end));
tab_crps_dmpv_each_full_ms = array2table(dmpv_crps_post_nowcast_ms,'RowNames',ModNames(2:end));
tab_crps_dmpv_each_full_pca = array2table(dmpv_crps_post_nowcast_pca,'RowNames',ModNames(2:end));
tab_crps_dmpv_each_full_ucomb = array2table(dmpv_crps_post_nowcast_ucomb,'RowNames',ModNames(2:end));


%%%% WQS Scores
dm_wqs_pre_nowcast_ms= [];
dm_wqs_pre_nowcast_pca= [];
dm_wqs_pre_nowcast_ucomb= [];
dmpv_wqs_pre_nowcast_ms = [];
dmpv_wqs_pre_nowcast_pca = [];
dmpv_wqs_pre_nowcast_ucomb = [];
dm_wqs_post_nowcast_ms = [];
dm_wqs_post_nowcast_pca = [];
dm_wqs_post_nowcast_ucomb = [];
dmpv_wqs_post_nowcast_ms = [];
dmpv_wqs_post_nowcast_pca = [];
dmpv_wqs_post_nowcast_ucomb = [];
wqs_each_pre_ms=nan(nummod-1,nper);
wqs_each_pre_pca=nan(nummod-1,nper);
wqs_each_pre_ucomb=nan(nummod-1,nper);
wqs_each_post_ms=nan(nummod-1,nper);
wqs_each_post_pca=nan(nummod-1,nper);
wqs_each_post_ucomb=nan(nummod-1,nper);
for nn=2:numel(fm)
    dm_wqs_pre_nowcast_ms = [dm_wqs_pre_nowcast_ms;dm_wqs_pre_each_ms.(fm{nn})'];
    dm_wqs_pre_nowcast_pca = [dm_wqs_pre_nowcast_pca;dm_wqs_pre_each_pca.(fm{nn})'];
    dm_wqs_pre_nowcast_ucomb = [dm_wqs_pre_nowcast_ucomb;dm_wqs_pre_each_ucomb.(fm{nn})'];
    dmpv_wqs_pre_nowcast_ms = [dmpv_wqs_pre_nowcast_ms;dmpv_wqs_pre_each_ms.(fm{nn})'];
    dmpv_wqs_pre_nowcast_pca = [dmpv_wqs_pre_nowcast_pca;dmpv_wqs_pre_each_pca.(fm{nn})'];
    dmpv_wqs_pre_nowcast_ucomb = [dmpv_wqs_pre_nowcast_ucomb;dmpv_wqs_pre_each_ucomb.(fm{nn})'];
    dm_wqs_post_nowcast_ms = [dm_wqs_post_nowcast_ms;dm_wqs_post_each_ms.(fm{nn})'];
    dm_wqs_post_nowcast_pca = [dm_wqs_post_nowcast_pca;dm_wqs_post_each_pca.(fm{nn})'];
    dm_wqs_post_nowcast_ucomb = [dm_wqs_post_nowcast_ucomb;dm_wqs_post_each_ucomb.(fm{nn})'];
    dmpv_wqs_post_nowcast_ms = [dmpv_wqs_post_nowcast_ms;dmpv_wqs_post_each_ms.(fm{nn})'];
    dmpv_wqs_post_nowcast_pca = [dmpv_wqs_post_nowcast_pca;dmpv_wqs_post_each_pca.(fm{nn})'];
    dmpv_wqs_post_nowcast_ucomb= [dmpv_wqs_post_nowcast_ucomb;dmpv_wqs_post_each_ucomb.(fm{nn})'];
    for np=1:nper
    wqs_each_pre_ms(nn-1,np) = wqs_pre.(fm{nn})(np)/wqs_pre.(fm{1})(np);
    wqs_each_pre_pca(nn-1,np) = wqs_pre.(fm{nn})(np)/wqs_pre.(fm{length(fm)})(np);
    wqs_each_pre_ucomb(nn-1,np) = wqs_pre.(fm{nn})(np)/wqs_pre.(fm{17})(np);
    wqs_each_full_ms(nn-1,np) = wqs_full.(fm{nn})(np)/wqs_full.(fm{1})(np);
    wqs_each_full_pca(nn-1,np) = wqs_full.(fm{nn})(np)/wqs_full.(fm{length(fm)})(np);
    wqs_each_full_ucomb(nn-1,np) = wqs_full.(fm{nn})(np)/wqs_full.(fm{17})(np);
    end
end
tab_wqs_each_pre_ms = array2table(wqs_each_pre_ms,'RowNames',ModNames(2:end));
tab_wqs_each_pre_pca = array2table(wqs_each_pre_pca,'RowNames',ModNames(2:end));
tab_wqs_each_pre_ucomb = array2table(wqs_each_pre_ucomb,'RowNames',ModNames(2:end));
tab_wqs_each_full_ms = array2table(wqs_each_full_ms,'RowNames',ModNames(2:end));
tab_wqs_each_full_pca = array2table(wqs_each_full_pca,'RowNames',ModNames(2:end));
tab_wqs_each_full_ucomb = array2table(wqs_each_full_ucomb,'RowNames',ModNames(2:end));
tab_wqs_dm_each_pre_ms = array2table(dm_wqs_pre_nowcast_ms,'RowNames',ModNames(2:end));
tab_wqs_dm_each_pre_pca = array2table(dm_wqs_pre_nowcast_pca,'RowNames',ModNames(2:end));
tab_wqs_dm_each_pre_ucomb= array2table(dm_wqs_pre_nowcast_ucomb,'RowNames',ModNames(2:end));
tab_wqs_dm_each_full_ms = array2table(dm_wqs_post_nowcast_ms,'RowNames',ModNames(2:end));
tab_wqs_dm_each_full_pca = array2table(dm_wqs_post_nowcast_pca,'RowNames',ModNames(2:end));
tab_wqs_dm_each_full_ucomb = array2table(dm_wqs_post_nowcast_ucomb,'RowNames',ModNames(2:end));
tab_wqs_dmpv_each_pre_ms = array2table(dmpv_wqs_pre_nowcast_ms,'RowNames',ModNames(2:end));
tab_wqs_dmpv_each_pre_pca = array2table(dmpv_wqs_pre_nowcast_pca,'RowNames',ModNames(2:end));
tab_wqs_dmpv_each_pre_ucomb = array2table(dmpv_wqs_pre_nowcast_ucomb,'RowNames',ModNames(2:end));
tab_wqs_dmpv_each_full_ms = array2table(dmpv_wqs_post_nowcast_ms,'RowNames',ModNames(2:end));
tab_wqs_dmpv_each_full_pca = array2table(dmpv_wqs_post_nowcast_pca,'RowNames',ModNames(2:end));
tab_wqs_dmpv_each_full_ucomb = array2table(dmpv_wqs_post_nowcast_ucomb,'RowNames',ModNames(2:end));


%% Table 4: Scores Averaged across Nowcast Periods: Against T-SV-GIGG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================== Retrieve RMSE and CRPS + compute DM =================== %
fm = fieldnames(res);

%RMSE and CRPS, for each nowcast period
for nn=1:numel(fm)
%%%%%%%%%% RMSE and CRPS, average over all nowcast periods
    rmsfe_gfc.(fm{nn}) = std(res.(fm{nn})(:,nfor_gfc),0,2);
    rmsfe_tranq.(fm{nn}) = std(res.(fm{nn})(:,nfor_tranq),0,2);
    rmsfe_covid.(fm{nn}) = std(res.(fm{nn})(:,nfor_covid),0,2);
    rmsfe_pre.(fm{nn}) = std(res.(fm{nn})(:,1:52),0,2);
    rmsfe_full.(fm{nn}) = std(res.(fm{nn})(:,1:nfor_full),0,2);
    crps_gfc.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_gfc),2);
    crps_tranq.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_tranq),2);
    crps_covid.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_covid),2);
    crps_pre.(fm{nn}) = mean(crps.(fm{nn})(:,1:52),2);
    crps_full.(fm{nn}) = mean(crps.(fm{nn})(:,1:nfor_full),2);
    wqs_gfc.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_gfc),2);
    wqs_tranq.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_tranq),2);
    wqs_covid.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_covid),2);
    wqs_full.(fm{nn}) = mean(wqs.(fm{nn})(:,1:nfor_full),2);
    wqs_pre.(fm{nn}) = mean(wqs.(fm{nn})(:,1:52),2);
end

%%%%%%%%%% DM stats Point forecasts
for nn=[1:2 4:numel(fm)]
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_gfc_joint.(fm{nn}),dmpv_gfc_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{3})(:,nfor_gfc),[],1),reshape(res.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_tranq_joint.(fm{nn}),dmpv_tranq_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{3})(:,nfor_tranq),[],1),reshape(res.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_covid_joint.(fm{nn}),dmpv_covid_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{3})(:,nfor_covid),[],1),reshape(res.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_pre_joint.(fm{nn}),dmpv_pre_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{3})(:,1:52),[],1),reshape(res.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_full_joint.(fm{nn}),dmpv_full_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{3})(:,1:nfor_full),[],1),reshape(res.(fm{nn})(:,1:nfor_full),[],1));
end

%%%%%%%%%% DM stats Density forecasts
for nn=[1:2 4:numel(fm)]
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_gfc_joint.(fm{nn}),dmpv_crps_gfc_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{3})(:,nfor_gfc),[],1),reshape(crps.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_tranq_joint.(fm{nn}),dmpv_crps_tranq_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{3})(:,nfor_tranq),[],1),reshape(crps.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_covid_joint.(fm{nn}),dmpv_crps_covid_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{3})(:,nfor_covid),[],1),reshape(crps.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_pre_joint.(fm{nn}),dmpv_crps_pre_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{3})(:,1:52),[],1),reshape(crps.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_full_joint.(fm{nn}),dmpv_crps_full_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{3})(:,1:nfor_full),[],1),reshape(crps.(fm{nn})(:,1:nfor_full),[],1));
end

%%%%%%%%%% DM stats WQS forecasts
for nn=[1:2 4:numel(fm)]
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_gfc_joint.(fm{nn}),dmpv_wqs_gfc_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{3})(:,nfor_gfc),[],1),reshape(wqs.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_tranq_joint.(fm{nn}),dmpv_wqs_tranq_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{3})(:,nfor_tranq),[],1),reshape(wqs.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_covid_joint.(fm{nn}),dmpv_wqs_covid_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{3})(:,nfor_covid),[],1),reshape(wqs.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_pre_joint.(fm{nn}),dmpv_wqs_pre_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{3})(:,1:52),[],1),reshape(wqs.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_full_joint.(fm{nn}),dmpv_wqs_full_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{3})(:,1:nfor_full),[],1),reshape(wqs.(fm{nn})(:,1:nfor_full),[],1));
end

temp = [];
for nn=[1:2 4:numel(fm)]
        temp =[temp; [mean(rmsfe_gfc.(fm{nn}))/mean(rmsfe_gfc.(fm{3})),mean(rmsfe_tranq.(fm{nn}))/mean(rmsfe_tranq.(fm{3})),mean(rmsfe_covid.(fm{nn}))/mean(rmsfe_covid.(fm{3})), mean(rmsfe_pre.(fm{nn}))/mean(rmsfe_pre.(fm{3})), mean(rmsfe_full.(fm{nn}))/mean(rmsfe_full.(fm{3})),...
        mean(crps_gfc.(fm{nn}))/mean(crps_gfc.(fm{3})), mean(crps_tranq.(fm{nn}))/mean(crps_tranq.(fm{3})),mean(crps_covid.(fm{nn}))/mean(crps_covid.(fm{3})),mean(crps_pre.(fm{nn}))/mean(crps_pre.(fm{3})),mean(crps_full.(fm{nn}))/mean(crps_full.(fm{3})),...
        mean(wqs_gfc.(fm{nn}))/mean(wqs_gfc.(fm{3})), mean(wqs_tranq.(fm{nn}))/mean(wqs_tranq.(fm{3})),mean(wqs_covid.(fm{nn}))/mean(wqs_covid.(fm{3})),mean(wqs_pre.(fm{nn}))/mean(wqs_pre.(fm{3})),mean(wqs_full.(fm{nn}))/mean(wqs_full.(fm{3}))],...
        dm_gfc_joint.(fm{nn}), dmpv_gfc_joint.(fm{nn}), ...
        dm_tranq_joint.(fm{nn}), dmpv_tranq_joint.(fm{nn}), ...
        dm_covid_joint.(fm{nn}), dmpv_covid_joint.(fm{nn}), ...
        dm_pre_joint.(fm{nn}), dmpv_pre_joint.(fm{nn}), ...
        dm_full_joint.(fm{nn}), dmpv_full_joint.(fm{nn}), ...
        dm_crps_gfc_joint.(fm{nn}), dmpv_crps_gfc_joint.(fm{nn}),...
        dm_crps_tranq_joint.(fm{nn}), dmpv_crps_tranq_joint.(fm{nn}),...
        dm_crps_covid_joint.(fm{nn}), dmpv_crps_covid_joint.(fm{nn}),...
        dm_crps_pre_joint.(fm{nn}), dmpv_crps_pre_joint.(fm{nn}),...
        dm_crps_full_joint.(fm{nn}), dmpv_crps_full_joint.(fm{nn}),...
        dm_wqs_gfc_joint.(fm{nn}), dmpv_wqs_gfc_joint.(fm{nn}),...
        dm_wqs_tranq_joint.(fm{nn}), dmpv_wqs_tranq_joint.(fm{nn}),...
        dm_wqs_covid_joint.(fm{nn}), dmpv_wqs_covid_joint.(fm{nn}),...
        dm_wqs_pre_joint.(fm{nn}), dmpv_wqs_pre_joint.(fm{nn}),...
        dm_wqs_full_joint.(fm{nn}), dmpv_wqs_full_joint.(fm{nn})];
end
%restable = [temp;[ mean(rmsfe_gfc.(fm{3})), mean(rmsfe_tranq.(fm{3})),mean(rmsfe_covid.(fm{3})),mean(rmsfe_pre.(fm{3})),mean(rmsfe_pre.(fm{3})),mean(crps_gfc.(fm{3})), mean(crps_tranq.(fm{3})),mean(crps_covid.(fm{3})),mean(crps_pre.(fm{3})),mean(crps_full.(fm{3})),nan(1,20)]];
restable = [temp(1:2,:);[ mean(rmsfe_gfc.(fm{3})), mean(rmsfe_tranq.(fm{3})),mean(rmsfe_covid.(fm{3})),mean(rmsfe_pre.(fm{3})),mean(rmsfe_pre.(fm{3})),mean(crps_gfc.(fm{3})), mean(crps_tranq.(fm{3})),mean(crps_covid.(fm{3})),mean(crps_pre.(fm{3})),mean(crps_full.(fm{3})),mean(wqs_tranq.(fm{3})),mean(wqs_covid.(fm{3})),mean(wqs_pre.(fm{3})),mean(wqs_full.(fm{3})),nan(1,31)];temp(3:end,:)];

Scores1 = ["RMSFE_GFC";"RMSFE_Tranquil";"RMSFE_Covid";"RMSFE_Pre";"RMSFE_Full";"CRPS_GFC";"CRPS_Tranquil";"CRPS_Covid";"CRPS_Pre";"CRPS_Full";"WQS_GFC";"WQS_Tranquil";"WQS_Covid";"WQS_Pre";"WQS_Full";"DM_joint_GFC";"DM_pv_joint_GFC";"DM_joint_Tranquil";"DM_pv_joint_Tranquil";"DM_joint_Covid";"DM_pv_joint_Covid";"DM_joint_Pre";"DM_pv_joint_Pre";"DM_joint_Full";"DM_pv_joint_Full";"DM_CRPS_joint_GFC";"DM_CRPS_pv_joint_GFC";"DM_CRPS_joint_Tranquil";"DM_CRPS_pv_joint_Tranquil";"DM_CRPS_joint_Covid";"DM_CRPS_pv_joint_Covid";"DM_CRPS_joint_Pre";"DM_CRPS_pv_joint_Pre";"DM_CRPS_joint_Full";"DM_CRPS_pv_joint_Full";"DM_WQS_joint_GFC";"DM_WQS_pv_joint_GFC";"DM_WQS_joint_Tranquil";"DM_WQS_pv_joint_Tranquil";"DM_WQS_joint_Covid";"DM_WQS_pv_joint_Covid";"DM_WQS_joint_Pre";"DM_WQS_pv_joint_Pre";"DM_WQS_joint_Full";"DM_WQS_pv_joint_Full"];
tab4 = array2table(restable,'VariableNames',Scores1,'RowNames',ModNames);


%% Table 4: Scores Averaged across Nowcast Periods: Against U-MIDAS-Comb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================== Retrieve RMSE and CRPS + compute DM =================== %
fm = fieldnames(res);

%RMSE and CRPS, for each nowcast period
for nn=[1:16 18:numel(fm)]
%%%%%%%%%% RMSE and CRPS, average over all nowcast periods
    rmsfe_gfc.(fm{nn}) = std(res.(fm{nn})(:,nfor_gfc),0,2);
    rmsfe_tranq.(fm{nn}) = std(res.(fm{nn})(:,nfor_tranq),0,2);
    rmsfe_covid.(fm{nn}) = std(res.(fm{nn})(:,nfor_covid),0,2);
    rmsfe_pre.(fm{nn}) = std(res.(fm{nn})(:,1:52),0,2);
    rmsfe_full.(fm{nn}) = std(res.(fm{nn})(:,1:nfor_full),0,2);
    crps_gfc.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_gfc),2);
    crps_tranq.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_tranq),2);
    crps_covid.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_covid),2);
    crps_pre.(fm{nn}) = mean(crps.(fm{nn})(:,1:52),2);
    crps_full.(fm{nn}) = mean(crps.(fm{nn})(:,1:nfor_full),2);
    wqs_gfc.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_gfc),2);
    wqs_tranq.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_tranq),2);
    wqs_covid.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_covid),2);
    wqs_pre.(fm{nn}) = mean(wqs.(fm{nn})(:,1:52),2);
    wqs_full.(fm{nn}) = mean(wqs.(fm{nn})(:,1:nfor_full),2);
end

%%%%%%%%%% DM stats Point forecasts
for nn=[1:16 18:numel(fm)]
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_gfc_joint.(fm{nn}),dmpv_gfc_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{17})(:,nfor_gfc),[],1),reshape(res.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_tranq_joint.(fm{nn}),dmpv_tranq_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{17})(:,nfor_tranq),[],1),reshape(res.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_covid_joint.(fm{nn}),dmpv_covid_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{17})(:,nfor_covid),[],1),reshape(res.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_pre_joint.(fm{nn}),dmpv_pre_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{17})(:,1:52),[],1),reshape(res.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_full_joint.(fm{nn}),dmpv_full_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{17})(:,1:nfor_full),[],1),reshape(res.(fm{nn})(:,1:nfor_full),[],1));
end

%%%%%%%%%% DM stats Density forecasts
for nn=[1:16 18:numel(fm)]
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_gfc_joint.(fm{nn}),dmpv_crps_gfc_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{17})(:,nfor_gfc),[],1),reshape(crps.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_tranq_joint.(fm{nn}),dmpv_crps_tranq_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{17})(:,nfor_tranq),[],1),reshape(crps.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_covid_joint.(fm{nn}),dmpv_crps_covid_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{17})(:,nfor_covid),[],1),reshape(crps.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_pre_joint.(fm{nn}),dmpv_crps_pre_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{17})(:,1:52),[],1),reshape(crps.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_full_joint.(fm{nn}),dmpv_crps_full_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{17})(:,1:nfor_full),[],1),reshape(crps.(fm{nn})(:,1:nfor_full),[],1));
end

%%%%%%%%%% DM stats Density forecasts
for nn=[1:16 18:numel(fm)]
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_gfc_joint.(fm{nn}),dmpv_wqs_gfc_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{17})(:,nfor_gfc),[],1),reshape(wqs.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_tranq_joint.(fm{nn}),dmpv_wqs_tranq_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{17})(:,nfor_tranq),[],1),reshape(wqs.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_covid_joint.(fm{nn}),dmpv_wqs_covid_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{17})(:,nfor_covid),[],1),reshape(wqs.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_pre_joint.(fm{nn}),dmpv_wqs_pre_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{17})(:,1:52),[],1),reshape(wqs.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_full_joint.(fm{nn}),dmpv_wqs_full_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{17})(:,1:nfor_full),[],1),reshape(wqs.(fm{nn})(:,1:nfor_full),[],1));
end

temp = [];
for nn=[1:16 18:numel(fm)]
        temp =[temp; [mean(rmsfe_gfc.(fm{nn}))/mean(rmsfe_gfc.(fm{17})),mean(rmsfe_tranq.(fm{nn}))/mean(rmsfe_tranq.(fm{17})),mean(rmsfe_covid.(fm{nn}))/mean(rmsfe_covid.(fm{17})), mean(rmsfe_pre.(fm{nn}))/mean(rmsfe_pre.(fm{17})), mean(rmsfe_full.(fm{nn}))/mean(rmsfe_full.(fm{17})),...
        mean(crps_gfc.(fm{nn}))/mean(crps_gfc.(fm{17})), mean(crps_tranq.(fm{nn}))/mean(crps_tranq.(fm{17})),mean(crps_covid.(fm{nn}))/mean(crps_covid.(fm{17})),mean(crps_pre.(fm{nn}))/mean(crps_pre.(fm{17})),mean(crps_full.(fm{nn}))/mean(crps_full.(fm{17})),...
        mean(wqs_gfc.(fm{nn}))/mean(wqs_gfc.(fm{17})), mean(wqs_tranq.(fm{nn}))/mean(wqs_tranq.(fm{17})),mean(wqs_covid.(fm{nn}))/mean(wqs_covid.(fm{17})),mean(wqs_pre.(fm{nn}))/mean(wqs_pre.(fm{17})),mean(wqs_full.(fm{nn}))/mean(wqs_full.(fm{17}))],...
        dm_gfc_joint.(fm{nn}), dmpv_gfc_joint.(fm{nn}), ...
        dm_tranq_joint.(fm{nn}), dmpv_tranq_joint.(fm{nn}), ...
        dm_covid_joint.(fm{nn}), dmpv_covid_joint.(fm{nn}), ...
        dm_pre_joint.(fm{nn}), dmpv_pre_joint.(fm{nn}), ...
        dm_full_joint.(fm{nn}), dmpv_full_joint.(fm{nn}), ...
        dm_crps_gfc_joint.(fm{nn}), dmpv_crps_gfc_joint.(fm{nn}),...
        dm_crps_tranq_joint.(fm{nn}), dmpv_crps_tranq_joint.(fm{nn}),...
        dm_crps_covid_joint.(fm{nn}), dmpv_crps_covid_joint.(fm{nn}),...
        dm_crps_pre_joint.(fm{nn}), dmpv_crps_pre_joint.(fm{nn}),...
        dm_crps_full_joint.(fm{nn}), dmpv_crps_full_joint.(fm{nn}),...
        dm_wqs_gfc_joint.(fm{nn}), dmpv_wqs_gfc_joint.(fm{nn}),...
        dm_wqs_tranq_joint.(fm{nn}), dmpv_wqs_tranq_joint.(fm{nn}),...
        dm_wqs_covid_joint.(fm{nn}), dmpv_wqs_covid_joint.(fm{nn}),...
        dm_wqs_pre_joint.(fm{nn}), dmpv_wqs_pre_joint.(fm{nn}),...
        dm_wqs_full_joint.(fm{nn}), dmpv_wqs_full_joint.(fm{nn})];
end
%restable = [temp;[ mean(rmsfe_gfc.(fm{17})), mean(rmsfe_tranq.(fm{17})),mean(rmsfe_covid.(fm{17})),mean(rmsfe_pre.(fm{17})),mean(rmsfe_pre.(fm{17})),mean(crps_gfc.(fm{17})), mean(crps_tranq.(fm{17})),mean(crps_covid.(fm{17})),mean(crps_pre.(fm{17})),mean(crps_full.(fm{17})),nan(1,20)]];
restable = [temp(1:16,:);[ mean(rmsfe_gfc.(fm{17})), mean(rmsfe_tranq.(fm{17})),mean(rmsfe_covid.(fm{17})),mean(rmsfe_pre.(fm{17})),mean(rmsfe_pre.(fm{17})),mean(crps_gfc.(fm{17})), mean(crps_tranq.(fm{17})),mean(crps_covid.(fm{17})),mean(crps_pre.(fm{17})),mean(crps_full.(fm{17})),mean(wqs_gfc.(fm{17})), mean(wqs_tranq.(fm{17})),mean(wqs_covid.(fm{17})),mean(wqs_pre.(fm{17})),mean(wqs_full.(fm{17})),nan(1,30)];temp(17:end,:)];

Scores1 = ["RMSFE_GFC";"RMSFE_Tranquil";"RMSFE_Covid";"RMSFE_Pre";"RMSFE_Full";"CRPS_GFC";"CRPS_Tranquil";"CRPS_Covid";"CRPS_Pre";"CRPS_Full";"WQS_GFC";"WQS_Tranquil";"WQS_Covid";"WQS_Pre";"WQS_Full";"DM_joint_GFC";"DM_pv_joint_GFC";"DM_joint_Tranquil";"DM_pv_joint_Tranquil";"DM_joint_Covid";"DM_pv_joint_Covid";"DM_joint_Pre";"DM_pv_joint_Pre";"DM_joint_Full";"DM_pv_joint_Full";"DM_CRPS_joint_GFC";"DM_CRPS_pv_joint_GFC";"DM_CRPS_joint_Tranquil";"DM_CRPS_pv_joint_Tranquil";"DM_CRPS_joint_Covid";"DM_CRPS_pv_joint_Covid";"DM_CRPS_joint_Pre";"DM_CRPS_pv_joint_Pre";"DM_CRPS_joint_Full";"DM_CRPS_pv_joint_Full";"DM_WQS_joint_GFC";"DM_WQS_pv_joint_GFC";"DM_WQS_joint_Tranquil";"DM_WQS_pv_joint_Tranquil";"DM_WQS_joint_Covid";"DM_WQS_pv_joint_Covid";"DM_WQS_joint_Pre";"DM_WQS_pv_joint_Pre";"DM_WQS_joint_Full";"DM_WQS_pv_joint_Full"];
tab5 = array2table(restable,'VariableNames',Scores1,'RowNames',ModNames);


%% Table 2: Scores Averaged across Nowcast Periods: Against PCA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================== Retrieve RMSE and CRPS + compute DM =================== %
fm = fieldnames(res);

%RMSE and CRPS, for each nowcast period
for nn=1:numel(fm)
%%%%%%%%%% RMSE and CRPS, average over all nowcast periods
    rmsfe_gfc.(fm{nn}) = std(res.(fm{nn})(:,nfor_gfc),0,2);
    rmsfe_tranq.(fm{nn}) = std(res.(fm{nn})(:,nfor_tranq),0,2);
    rmsfe_covid.(fm{nn}) = std(res.(fm{nn})(:,nfor_covid),0,2);
    rmsfe_pre.(fm{nn}) = std(res.(fm{nn})(:,1:52),0,2);
    rmsfe_full.(fm{nn}) = std(res.(fm{nn})(:,1:nfor_full),0,2);
    crps_gfc.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_gfc),2);
    crps_tranq.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_tranq),2);
    crps_covid.(fm{nn}) = mean(crps.(fm{nn})(:,nfor_covid),2);
    crps_pre.(fm{nn}) = mean(crps.(fm{nn})(:,1:52),2);
    crps_full.(fm{nn}) = mean(crps.(fm{nn})(:,1:nfor_full),2);
    wqs_gfc.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_gfc),2);
    wqs_tranq.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_tranq),2);
    wqs_covid.(fm{nn}) = mean(wqs.(fm{nn})(:,nfor_covid),2);
    wqs_pre.(fm{nn}) = mean(wqs.(fm{nn})(:,1:52),2);
    wqs_full.(fm{nn}) = mean(wqs.(fm{nn})(:,1:nfor_full),2);
end

%%%%%%%%%% DM stats Point forecasts
for nn=1:numel(fm)-1
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_gfc_joint.(fm{nn}),dmpv_gfc_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{length(fm)})(:,nfor_gfc),[],1),reshape(res.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_tranq_joint.(fm{nn}),dmpv_tranq_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{length(fm)})(:,nfor_tranq),[],1),reshape(res.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_covid_joint.(fm{nn}),dmpv_covid_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{length(fm)})(:,nfor_covid),[],1),reshape(res.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_pre_joint.(fm{nn}),dmpv_pre_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{length(fm)})(:,1:52),[],1),reshape(res.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_full_joint.(fm{nn}),dmpv_full_joint.(fm{nn})]= dmtest_modified(reshape(res.(fm{length(fm)})(:,1:nfor_full),[],1),reshape(res.(fm{nn})(:,1:nfor_full),[],1));
end

%%%%%%%%%% DM stats Density forecasts
for nn=1:numel(fm)-1
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_gfc_joint.(fm{nn}),dmpv_crps_gfc_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{length(fm)})(:,nfor_gfc),[],1),reshape(crps.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_tranq_joint.(fm{nn}),dmpv_crps_tranq_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{length(fm)})(:,nfor_tranq),[],1),reshape(crps.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_covid_joint.(fm{nn}),dmpv_crps_covid_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{length(fm)})(:,nfor_covid),[],1),reshape(crps.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_pre_joint.(fm{nn}),dmpv_crps_pre_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{length(fm)})(:,1:52),[],1),reshape(crps.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_crps_full_joint.(fm{nn}),dmpv_crps_full_joint.(fm{nn})]= dmtest_modified(reshape(crps.(fm{length(fm)})(:,1:nfor_full),[],1),reshape(crps.(fm{nn})(:,1:nfor_full),[],1));
end

%%%%%%%%%% WQS stats Density forecasts
for nn=1:numel(fm)-1
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_gfc_joint.(fm{nn}),dmpv_wqs_gfc_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{length(fm)})(:,nfor_gfc),[],1),reshape(wqs.(fm{nn})(:,nfor_gfc),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_tranq_joint.(fm{nn}),dmpv_wqs_tranq_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{length(fm)})(:,nfor_tranq),[],1),reshape(wqs.(fm{nn})(:,nfor_tranq),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_covid_joint.(fm{nn}),dmpv_wqs_covid_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{length(fm)})(:,nfor_covid),[],1),reshape(wqs.(fm{nn})(:,nfor_covid),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_pre_joint.(fm{nn}),dmpv_wqs_pre_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{length(fm)})(:,1:52),[],1),reshape(wqs.(fm{nn})(:,1:52),[],1));
    %%% over all quarters pre-pandemic and nowcast periods
    [dm_wqs_full_joint.(fm{nn}),dmpv_wqs_full_joint.(fm{nn})]= dmtest_modified(reshape(wqs.(fm{length(fm)})(:,1:nfor_full),[],1),reshape(wqs.(fm{nn})(:,1:nfor_full),[],1));
end

temp = [];
for nn=1:numel(fm)-1
        temp =[temp; [mean(rmsfe_gfc.(fm{nn}))/mean(rmsfe_gfc.(fm{length(fm)})),mean(rmsfe_tranq.(fm{nn}))/mean(rmsfe_tranq.(fm{length(fm)})),mean(rmsfe_covid.(fm{nn}))/mean(rmsfe_covid.(fm{length(fm)})), mean(rmsfe_pre.(fm{nn}))/mean(rmsfe_pre.(fm{length(fm)})), mean(rmsfe_full.(fm{nn}))/mean(rmsfe_full.(fm{length(fm)})),...
        mean(crps_gfc.(fm{nn}))/mean(crps_gfc.(fm{length(fm)})), mean(crps_tranq.(fm{nn}))/mean(crps_tranq.(fm{length(fm)})),mean(crps_covid.(fm{nn}))/mean(crps_covid.(fm{length(fm)})),mean(crps_pre.(fm{nn}))/mean(crps_pre.(fm{length(fm)})),mean(crps_full.(fm{nn}))/mean(crps_full.(fm{length(fm)})),...
        mean(wqs_gfc.(fm{nn}))/mean(wqs_gfc.(fm{length(fm)})), mean(wqs_tranq.(fm{nn}))/mean(wqs_tranq.(fm{length(fm)})),mean(wqs_covid.(fm{nn}))/mean(wqs_covid.(fm{length(fm)})),mean(wqs_pre.(fm{nn}))/mean(wqs_pre.(fm{length(fm)})),mean(wqs_full.(fm{nn}))/mean(wqs_full.(fm{length(fm)}))],...
        dm_gfc_joint.(fm{nn}), dmpv_gfc_joint.(fm{nn}), ...
        dm_tranq_joint.(fm{nn}), dmpv_tranq_joint.(fm{nn}), ...
        dm_covid_joint.(fm{nn}), dmpv_covid_joint.(fm{nn}), ...
        dm_pre_joint.(fm{nn}), dmpv_pre_joint.(fm{nn}), ...
        dm_full_joint.(fm{nn}), dmpv_full_joint.(fm{nn}), ...
        dm_crps_gfc_joint.(fm{nn}), dmpv_crps_gfc_joint.(fm{nn}),...
        dm_crps_tranq_joint.(fm{nn}), dmpv_crps_tranq_joint.(fm{nn}),...
        dm_crps_covid_joint.(fm{nn}), dmpv_crps_covid_joint.(fm{nn}),...
        dm_crps_pre_joint.(fm{nn}), dmpv_crps_pre_joint.(fm{nn}),...
        dm_crps_full_joint.(fm{nn}), dmpv_crps_full_joint.(fm{nn}),...
        dm_wqs_gfc_joint.(fm{nn}), dmpv_wqs_gfc_joint.(fm{nn}),...
        dm_wqs_tranq_joint.(fm{nn}), dmpv_wqs_tranq_joint.(fm{nn}),...
        dm_wqs_covid_joint.(fm{nn}), dmpv_wqs_covid_joint.(fm{nn}),...
        dm_wqs_pre_joint.(fm{nn}), dmpv_wqs_pre_joint.(fm{nn}),...
        dm_wqs_full_joint.(fm{nn}), dmpv_wqs_full_joint.(fm{nn})];
end
restable = [temp;[ mean(rmsfe_gfc.(fm{length(fm)})), mean(rmsfe_tranq.(fm{length(fm)})),mean(rmsfe_covid.(fm{length(fm)})),mean(rmsfe_pre.(fm{length(fm)})),mean(rmsfe_pre.(fm{length(fm)})),mean(crps_gfc.(fm{length(fm)})), mean(crps_tranq.(fm{length(fm)})),mean(crps_covid.(fm{length(fm)})),mean(crps_pre.(fm{length(fm)})),mean(crps_full.(fm{length(fm)})),mean(wqs_gfc.(fm{length(fm)})), mean(wqs_tranq.(fm{length(fm)})),mean(wqs_covid.(fm{length(fm)})),mean(wqs_pre.(fm{length(fm)})),mean(wqs_full.(fm{length(fm)})),nan(1,30)]];

Scores1 = ["RMSFE_GFC";"RMSFE_Tranquil";"RMSFE_Covid";"RMSFE_Pre";"RMSFE_Full";"CRPS_GFC";"CRPS_Tranquil";"CRPS_Covid";"CRPS_Pre";"CRPS_Full";"WQS_GFC";"WQS_Tranquil";"WQS_Covid";"WQS_Pre";"WQS_Full";"DM_joint_GFC";"DM_pv_joint_GFC";"DM_joint_Tranquil";"DM_pv_joint_Tranquil";"DM_joint_Covid";"DM_pv_joint_Covid";"DM_joint_Pre";"DM_pv_joint_Pre";"DM_joint_Full";"DM_pv_joint_Full";"DM_CRPS_joint_GFC";"DM_CRPS_pv_joint_GFC";"DM_CRPS_joint_Tranquil";"DM_CRPS_pv_joint_Tranquil";"DM_CRPS_joint_Covid";"DM_CRPS_pv_joint_Covid";"DM_CRPS_joint_Pre";"DM_CRPS_pv_joint_Pre";"DM_CRPS_joint_Full";"DM_CRPS_pv_joint_Full";"DM_WQS_joint_GFC";"DM_WQS_pv_joint_GFC";"DM_WQS_joint_Tranquil";"DM_WQS_pv_joint_Tranquil";"DM_WQS_joint_Covid";"DM_WQS_pv_joint_Covid";"DM_WQS_joint_Pre";"DM_WQS_pv_joint_Pre";"DM_WQS_joint_Full";"DM_WQS_pv_joint_Full"];
tab2 = array2table(restable,'VariableNames',Scores1,'RowNames',ModNames);


%% Convert to Excel files
loc = 'Output_iterated/';
tabloc = strcat(loc,'dm_results_fullset_iteratedpredictions_',date,'.xlsx');
writetable(tab1,tabloc,'Sheet','Average Scores - vs MS ','WriteRowNames',true)
writetable(tab2,tabloc,'Sheet','Average Scores - vs PCA ','WriteRowNames',true)
writetable(tab4,tabloc,'Sheet','Average Scores - vs Trend-SV ','WriteRowNames',true)
writetable(tab5,tabloc,'Sheet','Average Scores - vs UComb ','WriteRowNames',true)

writetable(tab3,tabloc,'Sheet','Average Scores per Subsample','WriteRowNames',true)

writetable(tab_rmsfe_full,tabloc,'Sheet','RMSFE each Full','WriteRowNames',true)
writetable(tab_rmsfe_gfc,tabloc,'Sheet','RMSFE each GFC','WriteRowNames',true)
writetable(tab_rmsfe_tranq,tabloc,'Sheet','RMSFE each Tranquil','WriteRowNames',true)
writetable(tab_rmsfe_covid,tabloc,'Sheet','RMSFE each Covid','WriteRowNames',true)
writetable(tab_rmsfe_pre,tabloc,'Sheet','RMSFE each Pre','WriteRowNames',true)

writetable(tab_crps_full,tabloc,'Sheet','CRPS each Full','WriteRowNames',true)
writetable(tab_crps_gfc,tabloc,'Sheet','CRPS each GFC','WriteRowNames',true)
writetable(tab_crps_tranq,tabloc,'Sheet','CRPS each Tranquil','WriteRowNames',true)
writetable(tab_crps_covid,tabloc,'Sheet','CRPS each Covid','WriteRowNames',true)
writetable(tab_crps_pre,tabloc,'Sheet','CRPS each Pre','WriteRowNames',true)

writetable(tab_wqs_full,tabloc,'Sheet','wqs each Full','WriteRowNames',true)
writetable(tab_wqs_gfc,tabloc,'Sheet','wqs each GFC','WriteRowNames',true)
writetable(tab_wqs_tranq,tabloc,'Sheet','wqs each Tranquil','WriteRowNames',true)
writetable(tab_wqs_covid,tabloc,'Sheet','wqs each Covid','WriteRowNames',true)
writetable(tab_wqs_pre,tabloc,'Sheet','wqs each Pre','WriteRowNames',true)


writetable(tab_dmpv_each_pre_ms,tabloc,'Sheet','DMpv pre each MS','WriteRowNames',true)
writetable(tab_dm_each_pre_ms,tabloc,'Sheet','DM pre each MS','WriteRowNames',true)
writetable(tab_rmse_each_pre_ms,tabloc,'Sheet','RMSE pre each MS','WriteRowNames',true)
writetable(tab_dmpv_each_full_ms,tabloc,'Sheet','DMpv full each MS','WriteRowNames',true)
writetable(tab_dm_each_full_ms,tabloc,'Sheet','DM full each MS','WriteRowNames',true)
writetable(tab_rmse_each_full_ms,tabloc,'Sheet','RMSE full each MS','WriteRowNames',true)

writetable(tab_crps_dmpv_each_pre_ms,tabloc,'Sheet','CRPS DMpv pre each MS','WriteRowNames',true)
writetable(tab_crps_dm_each_pre_ms,tabloc,'Sheet','CRPS DM pre each MS','WriteRowNames',true)
writetable(tab_crps_each_pre_ms,tabloc,'Sheet','CRPS pre each MS','WriteRowNames',true)
writetable(tab_crps_dmpv_each_full_ms,tabloc,'Sheet','CRPS DMpv full each MS','WriteRowNames',true)
writetable(tab_crps_dm_each_full_ms,tabloc,'Sheet','CRPS DM full each MS','WriteRowNames',true)
writetable(tab_crps_each_full_ms,tabloc,'Sheet','CRPS full each MS','WriteRowNames',true)

writetable(tab_wqs_dmpv_each_pre_ms,tabloc,'Sheet','wqs DMpv pre each MS','WriteRowNames',true)
writetable(tab_wqs_dm_each_pre_ms,tabloc,'Sheet','wqs DM pre each MS','WriteRowNames',true)
writetable(tab_wqs_each_pre_ms,tabloc,'Sheet','wqs pre each MS','WriteRowNames',true)
writetable(tab_wqs_dmpv_each_full_ms,tabloc,'Sheet','wqs DMpv full each MS','WriteRowNames',true)
writetable(tab_wqs_dm_each_full_ms,tabloc,'Sheet','wqs DM full each MS','WriteRowNames',true)
writetable(tab_wqs_each_full_ms,tabloc,'Sheet','wqs full each MS','WriteRowNames',true)

writetable(tab_dmpv_each_pre_pca,tabloc,'Sheet','DMpv pre each PCA','WriteRowNames',true)
writetable(tab_dm_each_pre_pca,tabloc,'Sheet','DM pre each PCA','WriteRowNames',true)
writetable(tab_rmse_each_pre_pca,tabloc,'Sheet','RMSE pre each PCA','WriteRowNames',true)
writetable(tab_dmpv_each_full_pca,tabloc,'Sheet','DMpv full each PCA','WriteRowNames',true)
writetable(tab_dm_each_full_pca,tabloc,'Sheet','DM full each PCA','WriteRowNames',true)
writetable(tab_rmse_each_full_pca,tabloc,'Sheet','RMSE full each PCA','WriteRowNames',true)

writetable(tab_crps_dmpv_each_pre_pca,tabloc,'Sheet','CRPS DMpv pre each PCA','WriteRowNames',true)
writetable(tab_crps_dm_each_pre_pca,tabloc,'Sheet','CRPS DM pre each PCA','WriteRowNames',true)
writetable(tab_crps_each_pre_pca,tabloc,'Sheet','CRPS pre each PCA','WriteRowNames',true)
writetable(tab_crps_dmpv_each_full_pca,tabloc,'Sheet','CRPS DMpv full each PCA','WriteRowNames',true)
writetable(tab_crps_dm_each_full_pca,tabloc,'Sheet','CRPS DM full each PCA','WriteRowNames',true)
writetable(tab_crps_each_full_pca,tabloc,'Sheet','CRPS full each PCA','WriteRowNames',true)

writetable(tab_wqs_dmpv_each_pre_pca,tabloc,'Sheet','wqs DMpv pre each PCA','WriteRowNames',true)
writetable(tab_wqs_dm_each_pre_pca,tabloc,'Sheet','wqs DM pre each PCA','WriteRowNames',true)
writetable(tab_wqs_each_pre_pca,tabloc,'Sheet','wqs pre each PCA','WriteRowNames',true)
writetable(tab_wqs_dmpv_each_full_pca,tabloc,'Sheet','wqs DMpv full each PCA','WriteRowNames',true)
writetable(tab_wqs_dm_each_full_pca,tabloc,'Sheet','wqs DM full each PCA','WriteRowNames',true)
writetable(tab_wqs_each_full_pca,tabloc,'Sheet','wqs full each PCA','WriteRowNames',true)

writetable(tab_dmpv_each_pre_ucomb,tabloc,'Sheet','DMpv pre each UComb','WriteRowNames',true)
writetable(tab_dm_each_pre_ucomb,tabloc,'Sheet','DM pre each UComb','WriteRowNames',true)
writetable(tab_rmse_each_pre_ucomb,tabloc,'Sheet','RMSE pre each UComb','WriteRowNames',true)
writetable(tab_dmpv_each_full_ucomb,tabloc,'Sheet','DMpv full each UComb','WriteRowNames',true)
writetable(tab_dm_each_full_ucomb,tabloc,'Sheet','DM full each UComb','WriteRowNames',true)
writetable(tab_rmse_each_full_ucomb,tabloc,'Sheet','RMSE full each UComb','WriteRowNames',true)

writetable(tab_crps_dmpv_each_pre_ucomb,tabloc,'Sheet','CRPS DMpv pre each UComb','WriteRowNames',true)
writetable(tab_crps_dm_each_pre_ucomb,tabloc,'Sheet','CRPS DM pre each UComb','WriteRowNames',true)
writetable(tab_crps_each_pre_ucomb,tabloc,'Sheet','CRPS pre each UComb','WriteRowNames',true)
writetable(tab_crps_dmpv_each_full_ucomb,tabloc,'Sheet','CRPS DMpv full each UComb','WriteRowNames',true)
writetable(tab_crps_dm_each_full_ucomb,tabloc,'Sheet','CRPS DM full each UComb','WriteRowNames',true)
writetable(tab_crps_each_full_ucomb,tabloc,'Sheet','CRPS full each UComb','WriteRowNames',true)

writetable(tab_wqs_dmpv_each_pre_ucomb,tabloc,'Sheet','wqs DMpv pre each UComb','WriteRowNames',true)
writetable(tab_wqs_dm_each_pre_ucomb,tabloc,'Sheet','wqs DM pre each UComb','WriteRowNames',true)
writetable(tab_wqs_each_pre_ucomb,tabloc,'Sheet','wqs pre each UComb','WriteRowNames',true)
writetable(tab_wqs_dmpv_each_full_ucomb,tabloc,'Sheet','wqs DMpv full each UComb','WriteRowNames',true)
writetable(tab_wqs_dm_each_full_ucomb,tabloc,'Sheet','wqs DM full each UComb','WriteRowNames',true)
writetable(tab_wqs_each_full_ucomb,tabloc,'Sheet','wqs full each UComb','WriteRowNames',true)



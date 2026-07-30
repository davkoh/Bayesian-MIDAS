%% dm_eval_simplified.m  –  Compact forecast-evaluation script
%  ---------------------------------------------------------------
%  Produces a single Excel file with one summary sheet per benchmark
%  (score ratios + DM test statistics across subsamples).
%  ---------------------------------------------------------------

clear all
addpath('functions/')

%% ===================================================================
%%  1.  MODEL CONFIGURATION
%% ===================================================================
outputfolder = fullfile(cd, 'output');

%  Each row: { display_name,  model_folder_name,  'single'|'combination' }
% model_spec = { ...
%    'T-SV-GIGG-fixed',   'T_fixedSV_OBS_fixedSV_norm_BMIDAS_almon_gigg_fixed_postspars_yes',  'single'; ...
%    'T-SV-GIGG-hier-a',  'T_fixedSV_OBS_fixedSV_norm_BMIDAS_almon_gigg_hier_a_postspars_yes', 'single'; ...
%};

model_spec = { ...
    'T-SV-fixed-GIGG-hier-a',       'T_fixed_SV_OBS_fixed_SV_norm_BMIDAS_almon_gigg_hier_a_postspars_yes',  'single'; ...
    'T-SV-PC-GIGG-hier-a',       'T_PC_OBS_PC_norm_BMIDAS_almon_gigg_hier_a_postspars_yes',  'single'; ...
%%% Add more models here after running Main_new.m with different settings, e.g.:
    'HS',     'T_none_OBS_none_norm_BMIDAS_almon_horseshoe__postspars_yes',  'single'; ...
    'MS',          'T_none_OBS_none_norm_BMIDAS_almon_MAL__postspars_yes',         'single'; ...
    'GIGG-hier-a',          'T_none_OBS_none_norm_BMIDAS_almon_gigg_hier_a_postspars_yes',         'single'; ...
};

rows_to_drop = [];

benchmark_names = ["MS","GIGG-hier-a"];

%% ===================================================================
%%  2.  SUBSAMPLE DEFINITION (date-based)
%% ===================================================================
%  Use '' to auto-detect start/end from the evaluation period dates.

subsample_def = { ...
    'GFC',       '',           'Dec-2009'; ...
    'Tranquil',  'Jan-2010',   'Dec-2019'; ...
    'Covid',     'Jan-2020',   ''; ...
    'Pre',       '',           'Dec-2019'; ...
    'Full',      '',           ''; ...
};

subsample_names = string(subsample_def(:,1))';
nsub = size(subsample_def, 1);

% WQS settings
wqs_quantiles = 0.05:0.05:0.95;
wqs_weighting = 3;   % 1 = uniform, 2 = centre-weighted, 3 = tail-weighted

%% ===================================================================
%%  3.  LOAD DATA
%% ===================================================================
nummod = size(model_spec, 1);
ModNames = string(model_spec(:,1));

fprintf('Loading %d models ...\n', nummod);
dq_ref = [];
model_found = true(nummod, 1);

for mm = 1:nummod
    matfile = fullfile(outputfolder, model_spec{mm,2}, 'results.mat');
    if ~isfile(matfile)
        warning('Output file not found – skipping "%s":\n  %s', ModNames(mm), matfile);
        model_found(mm) = false;
        res{mm} = [];  crp{mm} = [];  wq{mm} = [];             %#ok<SAGROW>
        continue
    end

    raw = load(matfile);
    raw = raw.output;

    if isempty(dq_ref)
        dq_ref = raw.d_q;
        if isdatetime(dq_ref), dq_ref.Format = 'MMM-yyyy'; end
    end

    drop = rows_to_drop;
    if ~isempty(drop)
        drop(drop == Inf) = size(raw.resid_all, 1);
        raw.resid_all(drop, :) = [];
        raw.crps_all(drop, :)  = [];
    end

    res{mm} = raw.resid_all;                                    
    crp{mm} = raw.crps_all;                                     

    if strcmp(model_spec{mm,3}, 'combination')
        pred = squeeze(mean(raw.y_pred_all, 2));
    else
        pred = raw.y_pred_all;
    end
    wq{mm} = calculateWQS(pred, raw.yf, wqs_quantiles, wqs_weighting);  
    if ~isempty(drop), wq{mm}(drop, :) = []; end

    fprintf('  [%2d/%2d]  %s\n', mm, nummod, ModNames(mm));
end

if any(~model_found)
    fprintf('  Keeping %d of %d models (skipped %d missing).\n', ...
        sum(model_found), nummod, sum(~model_found));
    res = res(model_found);  crp = crp(model_found);  wq = wq(model_found);
    ModNames   = ModNames(model_found);
    model_spec = model_spec(model_found, :);
    nummod     = sum(model_found);
end

nper = size(res{1}, 1);

%% ===================================================================
%%  3b. RESOLVE SUBSAMPLES
%% ===================================================================
subsample_ranges = cell(nsub, 1);
for ss = 1:nsub
    start_str = subsample_def{ss, 2};
    end_str   = subsample_def{ss, 3};
    if isempty(start_str), dt_start = dq_ref(1);
    else, dt_start = datetime(start_str, 'InputFormat', 'MMM-yyyy'); end
    if isempty(end_str), dt_end = dq_ref(end);
    else, dt_end = dateshift(datetime(end_str, 'InputFormat', 'MMM-yyyy'), 'end', 'month'); end
    subsample_ranges{ss} = find(dq_ref >= dt_start & dq_ref <= dt_end);
    fprintf('Subsample "%s": %d quarters\n', subsample_def{ss,1}, numel(subsample_ranges{ss}));
end

%% ===================================================================
%%  4.  COMPUTE SCORES (averaged over nowcast periods)
%% ===================================================================
rmsfe_avg = nan(nummod, nsub);
crps_avg  = nan(nummod, nsub);
wqs_avg   = nan(nummod, nsub);

for mm = 1:nummod
    for ss = 1:nsub
        cols = subsample_ranges{ss};
        rmsfe_avg(mm, ss) = mean(std(res{mm}(:, cols), 0, 2));
        crps_avg(mm, ss)  = mean(mean(crp{mm}(:, cols), 2));
        wqs_avg(mm, ss)   = mean(mean(wq{mm}(:, cols), 2));
    end
end

%% ===================================================================
%%  5.  BUILD LONG-FORMAT TABLE (one per benchmark)
%% ===================================================================
%  Layout: rows = model × subsample, columns = RMSFE, CRPS, WQS,
%          DM_pv_RMSFE, DM_pv_CRPS, DM_pv_WQS
%  Benchmark rows show absolute levels; others show ratios vs benchmark.
%  Significance: * p<0.10, ** p<0.05, *** p<0.01 (ratio < 1 = improvement).

bench_idx = nan(1, numel(benchmark_names));
for bb = 1:numel(benchmark_names)
    idx = find(ModNames == benchmark_names(bb));
    if isempty(idx), error('Benchmark "%s" not found.', benchmark_names(bb)); end
    bench_idx(bb) = idx;
end

all_tabs = struct();

for bb = 1:numel(bench_idx)
    bi = bench_idx(bb);
    bname = benchmark_names(bb);
    fprintf('\n--- Summary vs benchmark: %s ---\n', bname);

    nrows = nummod * nsub;
    row_model = strings(nrows, 1);
    row_sub   = strings(nrows, 1);
    c_rmsfe   = strings(nrows, 1);
    c_crps    = strings(nrows, 1);
    c_wqs     = strings(nrows, 1);

    rr = 0;
    for mm = 1:nummod
        for ss = 1:nsub
            rr = rr + 1;
            row_model(rr) = ModNames(mm);
            row_sub(rr)   = subsample_names(ss);

            if mm == bi
                % Benchmark: absolute levels
                c_rmsfe(rr) = sprintf('%.4f', rmsfe_avg(mm, ss));
                c_crps(rr)  = sprintf('%.4f', crps_avg(mm, ss));
                c_wqs(rr)   = sprintf('%.4f', wqs_avg(mm, ss));
            else
                % Ratios
                r_rmsfe = rmsfe_avg(mm, ss) / rmsfe_avg(bi, ss);
                r_crps  = crps_avg(mm, ss)  / crps_avg(bi, ss);
                r_wqs   = wqs_avg(mm, ss)   / wqs_avg(bi, ss);

                % DM p-values
                cols = subsample_ranges{ss};
                [~, pv_r] = dmtest_modified(reshape(res{bi}(:,cols),[],1), reshape(res{mm}(:,cols),[],1));
                [~, pv_c] = dmtest_modified(reshape(crp{bi}(:,cols),[],1), reshape(crp{mm}(:,cols),[],1));
                [~, pv_w] = dmtest_modified(reshape(wq{bi}(:,cols),[],1),  reshape(wq{mm}(:,cols),[],1));

                c_rmsfe(rr) = sprintf('%.4f%s', r_rmsfe, sig_stars(r_rmsfe, pv_r));
                c_crps(rr)  = sprintf('%.4f%s', r_crps,  sig_stars(r_crps,  pv_c));
                c_wqs(rr)   = sprintf('%.4f%s', r_wqs,   sig_stars(r_wqs,   pv_w));
            end
        end
    end

    T = table(row_model, row_sub, c_rmsfe, c_crps, c_wqs, ...
              'VariableNames', {'Model','Subsample','RMSFE','CRPS','WQS'});

    tag = matlab.lang.makeValidName(bname);
    all_tabs.(sprintf('vs_%s', tag)) = T;
end

%% ===================================================================
%%  6.  EXPORT TO EXCEL
%% ===================================================================
tabloc = fullfile(outputfolder, sprintf('dm_results_%s.xlsx', datestr(now,'yyyy_mm_dd')));

sheet_names = fieldnames(all_tabs);
fprintf('\nWriting %d sheet(s) to:\n  %s\n', numel(sheet_names), tabloc);
for kk = 1:numel(sheet_names)
    sname = strrep(sheet_names{kk}, '_', ' ');
    if length(sname) > 31, sname = sname(1:31); end
    writetable(all_tabs.(sheet_names{kk}), tabloc, ...
               'Sheet', sname);
end

fprintf('Done.  (* p<0.10, ** p<0.05, *** p<0.01 for ratio<1)\n');

%% ===================================================================
%%  Helper: significance stars
%% ===================================================================
%function s = sig_stars(ratio, pv)
%    if ratio >= 1 || isnan(pv)
%        s = '';
%    elseif pv < 0.01
%        s = '***';
%    elseif pv < 0.05
%        s = '**';
%    elseif pv < 0.10
%        s = '*';
%    else
%        s = '';
%    end
%end

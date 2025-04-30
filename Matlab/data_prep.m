% This script takes all the data and mixed freqeuncy choices and produces the data for estimation

% Things that need to be included: 
    % 1. create a hash for the data choices
    % 2. load already created matlab data and compare the hash (if it exists)
    % 3. If hash is new, then re-create the data

%% Load data
[data_quarterly, names_q]= xlsread('UK_data_bmidas.xlsx','QuarterlyData','A4:e500');
[data_monthly, names_m]= xlsread('UK_data_bmidas.xlsx','MonthlyData','A4:aq2000');
%%% important: 
%  - excel sheet should be read with variable names 
%  - first data row for data_monthly are transformation indices (will be used in clean_data.m, line 27)

%% Some housekeeping
clearvars input % Ignore this
% TODO: make sure that the below does not get overwritten
input.mstart = -3; % Starting month for each nowcast cycle. E.g: choose -3 for start in March if the latest reference month of the quarter is June.
input.mend = 2; % Ending month for each nowcast cycle. E.g: choose 2 for ending nowcasting in August if the reference quarter is June.

%% Check if new data needs to be created?

% Generate a data hash to identify data choices
data_choices = struct();
data_choices.beg_s = beg_s;
data_choices.end_s = end_s;
data_choices.beg_eval_per = beg_eval_per;
data_choices.groups = struct('sur', sur, 'act', act, 'lab', lab, 'pr', pr, 'mon', mon, ...
                             'mort', mort, 'fin', fin, 'ie', ie, 'vis', vis);
% Create substructures to allow inputs with different dimensions
data_choices.variables = struct('Var', struct('value', Var, 'size', size(Var)), ...
                                'Varq', struct('value', Varq, 'size', size(Varq)), ...
                                'Vars', struct('value', Vars, 'size', size(Vars)), ...
                                'Vara', struct('value', Vara, 'size', size(Vara)), ...
                                'Varl', struct('value', Varl, 'size', size(Varl)), ...
                                'Varp', struct('value', Varp, 'size', size(Varp)), ...
                                'Varm', struct('value', Varm, 'size', size(Varm)), ...
                                'Varmt', struct('value', Varmt, 'size', size(Varmt)), ...
                                'Varf', struct('value', Varf, 'size', size(Varf)), ...
                                'Vari', struct('value', Vari, 'size', size(Vari)), ...
                                'Varv', struct('value', Varv, 'size', size(Varv)));
data_choices.transformations = struct('dyoy', dyoy, 'stand', stand);
data_choices.lag_structure = struct('mismatch', mismatch, 'monthvars', monthvars, ...
                                    'almonrest', almonrest, 'poly', poly);
data_choices.calendar = struct('Var_delay', Var_delay, 'Varq_delay', Varq_delay, ...
                               'Vars_delay', Vars_delay, 'Vara_delay', Vara_delay, ...
                               'Varl_delay', Varl_delay, 'Varp_delay', Varp_delay, ...
                               'Varm_delay', Varm_delay, 'Varmt_delay', Varmt_delay, ...
                               'Varf_delay', Varf_delay, 'Vari_delay', Vari_delay, ...
                               'Varv_delay', Varv_delay, 'Var_pubgroup', Var_pubgroup, ...
                               'Varq_pubgroup', Varq_pubgroup, 'Vars_pubgroup', Vars_pubgroup, ...
                               'Vara_pubgroup', Vara_pubgroup, 'Varl_pubgroup', Varl_pubgroup, ...
                               'Varp_pubgroup', Varp_pubgroup, 'Varm_pubgroup', Varm_pubgroup, ...
                               'Varmt_pubgroup', Varmt_pubgroup, 'Varf_pubgroup', Varf_pubgroup, ...
                               'Vari_pubgroup', Vari_pubgroup, 'Varv_pubgroup', Varv_pubgroup);


% Serialize the structure and compute a hash
data_choices_serialized = jsonencode(data_choices);
data_hash = DataHash(data_choices_serialized);

% Save the hash for comparison in the next step
hash_file = 'Data/data_hash.mat';
if isfile(hash_file)
    load(hash_file, 'previous_hash');
    if strcmp(data_hash, previous_hash)
        disp('Data choices unchanged. Skipping data preparation.');
    else
        disp('Data choices changed. Proceeding with data preparation.');
        save(hash_file, 'data_hash', '-v7.3');
    end
else
    disp('No previous hash found. Proceeding with data preparation.');
    save(hash_file, 'data_hash', '-v7.3');
end

% if the data exists, don't run the rest:
if strcmp(data_hash, previous_hash)
    load('UK_dat_2024.mat');
    disp('Loaded existing data.');
    return;
else 

%% Do the data handling
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

end
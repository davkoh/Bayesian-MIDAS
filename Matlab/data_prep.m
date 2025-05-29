% This script takes all the data and mixed freqeuncy choices and produces the data for estimation


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
% Check if variables exist before assigning them to the structure
if exist('sur', 'var'), data_choices.groups.sur = sur; else, data_choices.groups.sur = []; end
if exist('act', 'var'), data_choices.groups.act = act; else, data_choices.groups.act = []; end
if exist('lab', 'var'), data_choices.groups.lab = lab; else, data_choices.groups.lab = []; end
if exist('pr', 'var'), data_choices.groups.pr = pr; else, data_choices.groups.pr = []; end
if exist('mon', 'var'), data_choices.groups.mon = mon; else, data_choices.groups.mon = []; end
if exist('mort', 'var'), data_choices.groups.mort = mort; else, data_choices.groups.mort = []; end
if exist('fin', 'var'), data_choices.groups.fin = fin; else, data_choices.groups.fin = []; end
if exist('ie', 'var'), data_choices.groups.ie = ie; else, data_choices.groups.ie = []; end
if exist('vis', 'var'), data_choices.groups.vis = vis; else, data_choices.groups.vis = []; end

% Create substructures to allow inputs with different dimensions
data_choices.variables = struct();
variable_names = {'Var', 'Varq', 'Vars', 'Vara', 'Varl', 'Varp', 'Varm', 'Varmt', 'Varf', 'Vari', 'Varv'};
for i = 1:length(variable_names)
    var_name = variable_names{i};
    if exist(var_name, 'var')
        data_choices.variables.(var_name) = struct('value', eval(var_name), 'size', size(eval(var_name)));
    else
        data_choices.variables.(var_name) = struct('value', [], 'size', []);
    end
end
data_choices.transformations = struct('dyoy', dyoy, 'stand', stand);
data_choices.lag_structure = struct('mismatch', mismatch, 'monthvars', monthvars, ...
                                    'almonrest', almonrest, 'poly', poly);
calendar_fields = {'Var_delay', 'Varq_delay', 'Vars_delay', 'Vara_delay', 'Varl_delay', ...
                   'Varp_delay', 'Varm_delay', 'Varmt_delay', 'Varf_delay', 'Vari_delay', ...
                   'Varv_delay', 'Var_pubgroup', 'Varq_pubgroup', 'Vars_pubgroup', ...
                   'Vara_pubgroup', 'Varl_pubgroup', 'Varp_pubgroup', 'Varm_pubgroup', ...
                   'Varmt_pubgroup', 'Varf_pubgroup', 'Vari_pubgroup', 'Varv_pubgroup'};
data_choices.calendar = struct();
for i = 1:length(calendar_fields)
    field_name = calendar_fields{i};
    if exist(field_name, 'var')
        data_choices.calendar.(field_name) = eval(field_name);
    else
        data_choices.calendar.(field_name) = [];
    end
end


% Serialize the structure and compute a hash
data_choices_serialized = jsonencode(data_choices);
data_hash = DataHash(data_choices_serialized);

% Define the hash file path
hash_file = fullfile('Data', 'data_hash.mat');

% Check if the hash file exists
if isfile(hash_file)
    % Load the previous hash
    loaded_data = load(hash_file);
    if isfield(loaded_data, 'data_hash')
        previous_hash = loaded_data.data_hash;
        if strcmp(data_hash, previous_hash)
            disp('Data choices unchanged. Skipping data preparation.');
            % If the hash matches, load existing data and exit
            if isfile(fullfile('Data', 'UK_dat_2024.mat'))
                load(fullfile('Data', 'UK_dat_2024.mat'));
                disp('Loaded existing data.');
                return; % Skip further processing
            else
                disp('Data file not found. Proceeding with data preparation.');
            end
        else
            disp('Data choices changed. Proceeding with data preparation.');
        end
    else
        disp('Hash file not found. Proceeding with data preparation.');
    end
else
    disp('No previous hash found. Proceeding with data preparation.');
end

% Save the current hash for future comparison
save(hash_file, 'data_hash', '-v7.3');


[data_quarterly, names_q]= xlsread('MF_FAME_FULL.xlsx','QuarterlyData','A4:e500');
[data_monthly, names_m]= xlsread('MF_FAME_FULL.xlsx','MonthlyData','A4:aq2000');
%%% important: 
%  - excel sheet should be read with variable names 
%  - first data row for data_monthly are transformation indices (will be used in clean_data.m, line 27)

%%%%%%%% Does the data cleaning 
clean_data    %%%% transform and plot data, and prepare data for estimation  

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



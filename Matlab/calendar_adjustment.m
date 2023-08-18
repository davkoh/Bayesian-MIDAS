%% Publication Calendar Selection Accirding to Nowcast Script Choices

% This script loads in the pre-defined excel calendar, according to chosen
% indicators and amount of lags. Note that variables which are not included
% in the calendar cannot be included in the nowcast exercise.

% Publication calendar is organised into 1) nowcast period (starting row
% 5), 2) themes (group of indicators that get published together), 3)
% months in which the variables are published (start: 3 months before first
% month of reference quarter, end: 2 months after first month of reference
% quarter), 4) arrays which define the publication sequence for each
% monthly lag of the indicator (starting at last month of reference quarter
% denoted by 0 to 12 months back from that month denoted by 11, this then 
% encompasses 12 months worth of lags that can be possibly included). For
% any further questions, please contact david.kohns@aalto.fi

%% 1) Load in Publication Calendar
VarInfo = readtable('publication_calendar.xlsx','Sheet','Variables','Format','auto');
Varpub = readtable('publication_calendar.xlsx','Sheet','Calendar','Format','auto');

%% 2) According to Variable names and monthly lag length selection, take out correct publication schedules
idxtemp = find(strcmp('Name',Varpub.VariableID) == 1);
Names = Varpub{idxtemp,2:end};
Names =  Names(~cellfun('isempty',Names));
puball = Varpub{4:end,4:end};
puball = str2double(puball);

% Monthly indicators
pub_m = [];
for i = 1: size(y_m,2)
    pubtemp = puball(:,find( strcmp(varnames_qm{i},Names) == 1)); %varnames
    pubtemp = pubtemp(:,1:monthvars);
    pub_m = [pub_m pubtemp];
end

% Quarterly variables
if lags>0
pub_q = [];
for i = 1: size(varnames_qm,2)-size(y_m,2)
    pubtemp = puball(:,find( strcmp(varnames_qm{size(y_m,2)+i},Names) == 1)); %varnames
    pubtemp = pubtemp(:,1:lags);
    pub_q = [pub_q pubtemp];
end

puball = [pub_m pub_q];

else
    puball = pub_m;
    end

%% 3) Retrieve group structure for potential Almon restrictions
Varpub = readtable('publication_calendar.xlsx','Sheet','Calendar','ReadVariableNames',false,'Format','auto');
groupall = Varpub{1,4:end};
groupall = str2double(groupall);

group_m = [];
for i = 1:size(y_m,2)
    grouptemp = groupall(:,find( strcmp(varnames_qm{i},Names) == 1)); %varnames
    grouptemp = grouptemp(:,1:monthvars);
    group_m = [group_m  grouptemp];
end

if lags>0
group_q = [];
for i = 1: size(varnames_qm,2)-size(y_m,2)
    grouptemp = groupall(:,find( strcmp(varnames_qm{size(y_m,2)+i},Names) == 1)); %varnames
    grouptemp = grouptemp(:,1:lags);
    group_q = [group_q grouptemp];
end

groupall = [group_m group_q];

else
    groupall = group_m;
    end


clear Varpub idxtemp grouptemp VarInfo pubtemp

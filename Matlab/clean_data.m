

%% 1) select data series
hard_dummy = [];
if sur ==1,   Var = [Var; Vars]; hard_dummy = [hard_dummy,zeros(1,size(Vars,1))]; end
if act ==1,   Var = [Var; Vara]; hard_dummy = [hard_dummy,ones(1,size(Vara,1))]; end
if lab ==1,   Var = [Var; Varl]; hard_dummy = [hard_dummy,ones(1,size(Varl,1))];end
if pr  ==1,    Var = [Var; Varp]; hard_dummy = [hard_dummy,ones(1,size(Varp,1))]; end
if mon ==1,   Var = [Var; Varm]; hard_dummy = [hard_dummy,ones(1,size(Varm,1))]; end
if mort ==1,   Var = [Var; Varmt]; hard_dummy = [hard_dummy,ones(1,size(Varmt,1))]; end
if fin ==1,   Var = [Var; Varf]; hard_dummy = [hard_dummy,ones(1,size(Varf,1))]; end
if ie  ==1,    Var = [Var; Vari]; hard_dummy = [hard_dummy,zeros(1,size(Vari,1))]; end
if vis ==1,   Var = [Var; Varv]; hard_dummy = [hard_dummy,ones(1,size(Varv,1))]; end                                                          

% select correct series
Var = Var(~cellfun('isempty',Var));
data_incl_m = find_variable_indices(Var,strrep(names_m(1,2:end)',' ',''));
names_incl_m=names_m(1,1+data_incl_m);
Varq = Varq(~cellfun('isempty',Varq));
data_incl_q = find_variable_indices(Varq,strrep(names_q(1,2:end)',' ',''));
names_incl_q=names_q(1,1+data_incl_q);

%% 2) define data sets and names
data_m = data_monthly(2:end,data_incl_m);
transf = data_monthly(1,data_incl_m);        %transforms for monthly data defined directly in excel
[T,n] = size(data_m);
varnames = names_m(1:T+2,[1;data_incl_m+1]);
d_m = datetime(names_m(4:T+2,1)); % date names

data_q = data_quarterly(:,data_incl_q);
[Tq,n_q] = size(data_q);
varnames_q = names_q(1:Tq+2,[1;data_incl_q+1]);
d_q = datetime(names_q(3:Tq+2,1));

%%% variables included
if length(data_incl_q)==1   % quarterly series included
    varq = 'GDP';
elseif length(data_incl_q)==2
    varq = 'GDPCON';
elseif length(data_incl_q)>=3
    varq = 'GDPCONINVHOUR';
end
% elseif length(data_incl_q)==4
%     varq = 'GDPCONHOURINV';

varm = '';   % groups of monthly series included
if act==1 && sur==1 && act==1&&lab==1 && pr==1&& mon==1 && mort==1 && fin==1 && ie==1 && vis == 1
    varm = char([varm, '_allm']);
else
    if sur==1,  varm = char([varm,'_sur']); end
    if act==1,  varm = char([varm,'_act']); end
    if lab==1,  varm = char([varm,'_lab']); end
    if pr==1,   varm = char([varm,'_pr']); end
    if mon==1,  varm = char([varm,'_mon']); end
    if mort==1, varm = char([varm,'_mort']); end
    if fin==1,  varm = char([varm,'_fin']);end
    if ie==1,   varm = char([varm,'_ie']);end
    if vis==1,  varm = char([varm,'_vis']);end
end
if stand==1, varm = char([varm,'_std']);end


clear data_monthly data_quarterly names_m names_q Var Varq Varm Varp Varl Vara Vari Varv Varf Vars

%%  3)Transform data
data_m_tr = nan(size(data_m));
data_q_tr = nan(size(data_q));

for ii=1:n
    if transf(ii)==0      % no transform
        data_m_tr(:,ii)=data_m(:,ii)./sqrt(nanvar(data_m(:,ii)));
    elseif transf(ii)==1   % logs
        data_m_tr(:,ii)=log(data_m(:,ii)).*100;
    elseif transf(ii)==2   % first diff
        if dyoy==1
        data_m_tr(13:end,ii)=data_m(13:end,ii)-data_m(1:end-12,ii);
        else
        data_m_tr(2:end,ii)=diff(data_m(:,ii));
        end
    elseif transf(ii)==3   % log first diff
        if dyoy==1
            %data_m_tr(13:end,ii) = (log(data_m(13:Tm,:))-log(data_m(1:Tm-12,ii))).*100;
            data_m_tr(13:end,ii) =(data_m(13:T,ii)./data_m(1:T-12,ii))*100-100;         %%% or percentage change instead (can make a difference with large fluctuations)
        else
            %data_m_tr(2:end,ii)=diff(log(data_m(:,ii))).*100;
            data_m_tr(2:end,ii)= (data_m(2:T,ii)./data_m(1:T-1,ii))*100-100;             %%% or percentage change instead (can make a difference with large fluctuations)
        end
    end
    %%% standardise series
    if stand ==1 %|| transf(ii)==0 
        nanvars_m(ii) = nanvar(data_m_tr(:,ii)); %ADDED (RB) variance of monthly GDP used for rescaling
        if ii==1
    data_m_tr (:,ii) = data_m_tr(:,ii)./sqrt(nanvar(data_m_tr(:,ii))); %standardise series
        else
   data_m_tr (:,ii) = (data_m_tr(:,ii)-nanmean(data_m_tr(:,ii)))./sqrt(nanvar(data_m_tr(:,ii))); %standardise series
        end
    end
end

for ii=1:n_q
    if dyoy ==1    %y-o-y growth rates
        %data_q_tr(4:end,ii) = (log(data_q(5:Tm,:))-log(data_q(1:Tm-4,ii))).*100;
        data_q_tr(4:end,ii) = (data_q(5:Tq,ii)./data_q(1:Tq-4,ii))*100-100;
         %data_m_tr(13:end,ii)=(data_m(13:T,ii)./data_m(1:T-12,ii))*100-100;        %%% or percentage change instead (can make a difference with large fluctuations)
    else    % log first diff
        %data_q_tr(:,ii)=diff(log(data_q(:,ii))).*100;
         data_q_tr(2:end,ii)=(data_q(2:Tq,ii)./data_q(1:Tq-1,ii))*100-100;           %%% or percentage change instead (can make a difference with large fluctuations)
    end
    %%% standardise series
    if stand ==1
       % data_q_tr = data_q_tr./sqrt(nanvars_m(1)); %ADDED and next line uncommented (RB)
    data_m_tr = (data_m_tr-nanmean(data_m_tr))./sqrt(nanvar(data_m_tr));      %standardise series
    end
end

% drop NaNs at the beginning
if dyoy == 1
    data_m_tr = data_m_tr(13:end,:);
else
    data_m_tr = data_m_tr(2:end,:);
end


%% Select, Transform and plot data

if dyoy ==1 
    d_m = d_m(12:end,:);
end
t1_m = find(d_m==beg_s)-2;
t1_q = find(d_q==beg_s)+1;
tend_m = find(year(d_m)==year(end_s),1,'last');
tend_q = find(d_q==end_s);

y_q = data_q_tr(t1_q:tend_q,:);
y = data_m_tr(t1_m:tend_m,:);
T = length(y);
Tq = length(y_q);
d_m = d_m(t1_m:tend_m);
d_q = d_q(t1_q:tend_q);

y_m = y;
varnames_qm = [varnames(1,2:end),varnames_q(1,2:end)];

% Adjust awkward space in monthly hours worked series
if find(strcmp('Hours ',varnames_qm)==1) > 0
    idxhours = find(strcmp('Hours ',varnames_qm)==1);
varnames_qm{1,idxhours} = 'Hours';
end

clear make_plot data_incl_m data_incl_q data_m data_q data_qr_tr data_m_tr figj hard_dummy n n_q save_data save_res stand t1_m t1_q tend_m tend_q varm Varmt varnames varnames_1 varq y transf trans_q idxhours

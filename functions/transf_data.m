function [data_m_tr, data_q_tr] = transf_data(data_m, data_q, transf, dyoy, stand)

data_m_tr = nan(size(data_m));
data_q_tr = nan(size(data_q));
[T,n]=size(data_m);
[Tq,n_q]= size(data_q);

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

% drop NaNs at the beginning, but make sure we are rounded to quarters

% data_m_tr = data_m_tr(4+dyoy*9:end,:);
% data_q_tr= data_q_tr(2+dyoy*3:end,:);
data_m_tr = data_m_tr(2+dyoy*11:end,:);
data_q_tr= data_q_tr(1+dyoy*3:end,:);

end

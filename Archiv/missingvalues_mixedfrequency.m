%% Data Cleaning Helper Script for Nowcasting Script

% This sheets does 3 things: 1) Adjust monthly and quarterly data vectors
% to the user choices, 2) Fill up any missing values with a PCA approach 3)
% bring the monthly data vector into UMIDAS form. 

%% 1: Adjust monthly and quarter series according to user choices

% Adjust quarterly series for nowcasting sample
y_q0 = y_q; % Save instance of full length of quarterly variables fir lag creation
d_q0 = d_q;
y_q = y_q(find(dqstart == d_q):find(dqend == d_q),: );
d_q = d_q(find(dqstart == d_q):find(dqend == d_q) );

y = y_q(:,1);

% Adjust monthly series for nowcasting
y_m = y_m(find(mstart == d_m):find(mend == d_m),:); 
d_m = d_m(find(mstart == d_m):find(mend == d_m)); 

% Start: Get availability of latest data for nowcasting conditional on availability
avail_ind = ~isnan(y_m(find(mend == d_m)-monthvars+1:find(mend == d_m),:));
% End: Get availability of latest data for nowcasting conditional on availability


%% 2: Fill up missing values with PCA approach

% Obtain factors to interpolate missing monthly values
if size(y_m,2)<2
    y_m = fillmissing(y_m,'nearest');  
else
[~, PCAs , ~ , ~ , ~ , ~] = pca((y_m-nanmean(y_m))./(nanvar(y_m)),'algorithm','als'); % PCA with missing values


if size(y_m,2)<6
    f_m = PCAs(:,1:size(y_m,2)); % First PCA is Factor
else
f_m = PCAs(:,1:6); % First PCA is Factor
end
% Use factor to fill up variables which have missing values
NoNaN_m = isnan(y_m)==0; % Specify Vector that finds the missings in monthly

% First create factor
flag = lagmatrix(f_m,0:1-1);
i_nl = sum(isnan(flag),2)==0;
c = nanmean(y_m);
for i = 1:size(y_m,2)
    idx_i2 =(NoNaN_m(:,i).*i_nl)==1;
    F = flag(idx_i2,:);
    x_i = y_m(idx_i2==1,i);
    T_i = length(F);
    lam_m(:,i) = ([ones(size(F,1),1) F]'*[ones(size(F,1),1) F])\([ones(size(F,1),1) F]'*x_i); % loadings
end

% Now use factor to fill out missing values
y_mm = []; % storage for filled up monthly series
for i = 1:size(y_m,2)
id = i;
y_temp = y_m(:,i);
y_temp_fitted = [ones(size(f_m,1),1) f_m]*lam_m(:,id);
id_nan = isnan(y_m(1:end,id));
y_test = y_m(1:end,id);
y_test(find(id_nan==1),:) = y_temp_fitted(find(id_nan==1),:);
y_mm = [y_mm y_test];
end

y_m = y_mm;
end


%% 3: Transform to UMIDAS data 

Xm = NaN(size(y_q,1),monthvars*size(y_m,2));
idx = 1:monthvars:size(Xm,2);
mlags = monthvars/(mismatch)-1;
for j = 1:size(y_m,2)
    xtemp0 = [];
    for i = mismatch:-1:1
        xtemp0 = [xtemp0 y_m(i:mismatch:end,j)];
    end

    
    if mlags>0
    xtemp1 = [];
    xtemp1 = [xtemp1 xtemp0(1+mlags:end,:)];
    for i = 1:mlags
    xtemp1 = [xtemp1 xtemp0(1+mlags-i:end-i,:);];
    end
    else
        xtemp1 = xtemp0;
    end
    Xm(:,idx(j):idx(j)+monthvars-1)= xtemp1;
end

% Start: Skip-sample the availability indicator for nowcasting with most
% recent data

avail_ind_temp = nan(1,monthvars*size(y_m,2));

for i = 1:size(y_m,2)
    avail_ind_temp(i*monthvars-monthvars+1:i*monthvars) = avail_ind(i);
end

avail_ind = vec(avail_ind)';

% End: Skip-sample the availability indicator for nowcasting with most
% recent data

% Add quarterly variables

if lags>0
    idxlag = find( strcmp('GDP_Q',varnames_qm(size(y_m,2)+1:end)) == 1 );
    for i = 1:lags
        Xm = [Xm y_q0(find(dqstart == d_q0)-i:end-i,idxlag)];
    end
end

clear c avail_ind_temp data_q_tr F f_m flag i_nl id id_nan idx idx_i2 j lam_m NoNaN_m PCAs T_i transf_q varnames_q x_i xtemp0 xtemp1 y_mm y_q y_temp y_temp_fitted y_test y_q0 d_q0

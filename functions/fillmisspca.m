function [y_nomiss] = fillmisspca(y, numpca)
%%% compute Principal components 
%%% use first numpca components to fill missing values

% Obtain factors to interpolate missing monthly values
if size(y,2)<2
    y = fillmissing(y,'nearest');  
else
[~, PCAs , ~ , ~ , ~ , ~] = pca((y-nanmean(y))./(nanvar(y)),'algorithm','als'); % PCA with missing values

if size(y,2)<numpca
    f_m = PCAs(:,1:size(y,2)); % First PCA is Factor
else
f_m = PCAs(:,1:numpca); % Use only first 6 factors for interpolation
end
% Use factor to fill up variables which have missing values
NoNaN_m = ~isnan(y); % Specify Vector that finds the missings in monthly

% First create factor
flag = lagmatrix(f_m,0:1-1);
i_nl = sum(isnan(flag),2)==0;
c = nanmean(y);
for i = 1:size(y,2)
    idx_i2 =(NoNaN_m(:,i).*i_nl)==1;
    F = flag(idx_i2,:);
    x_i = y(idx_i2==1,i);
    T_i = length(F);
    lam_m(:,i) = ([ones(size(F,1),1) F]'*[ones(size(F,1),1) F])\([ones(size(F,1),1) F]'*x_i); % loadings
end

% Now use factor to fill out missing values
y_nm = []; % storage for filled up monthly series
for i = 1:size(y,2)
id = i;
y_temp = y(:,i);
y_temp_fitted = [ones(size(f_m,1),1) f_m]*lam_m(:,id);
id_nan = isnan(y(1:end,id));
y_test = y(1:end,id);
y_test(find(id_nan==1),:) = y_temp_fitted(find(id_nan==1),:);
y_nm = [y_nm y_test];
end

y_nomiss = y_nm;
end
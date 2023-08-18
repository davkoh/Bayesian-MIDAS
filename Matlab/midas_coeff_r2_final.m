function [beta] = midas_coeff_r2_final(a,X,grouping,polyt)
grall = unique(grouping); % unique groups due to lags of individual monthly variables.
r = 2; % number of restrictions, for this function, 2.
poly = polyt - r; % polynomial degree adjusted for restrictions
lagl = histc(grouping, grall); % lag length of each covariate;
idxt = [];
for i = 1:size(grall,2)
    temp = repmat(i,1,poly);
    idxt = [idxt temp];
end

beta = zeros( size(a,1),max(lagl),size(grall,2)); % Storage matrix for the tranformed coefficients

for m = 1:size(grall,2)
    idx = find(idxt == m);
    bt = a(:,idx);
    betatemp = zeros(size(a,1),lagl(m),poly);
for j = 1: poly
for i = 0:lagl(m)-1
    %betatemp(:,i,j) = bt(:,j)*(i^(j+1) + j*lagl(m)^(j+1) - (j+1)*lagl(m)^j*i); 
    %beta(i) = bt(:,j)*(i^(j+1) + j*lagl(m)^(j+1) - (j+1)*lagl(m)^j*i); 
    betatemp(:,i+1,j) = bt(:,j)*(i^(j+1) + j*(lagl(m)-1)^(j+1) - (j+1)*(lagl(m)-1)^j*(i)); 
end
end
beta(:,1:lagl(m),m) = squeeze(sum(betatemp,3));

end
%% Group-SAVS Function
function [pSAVS] = group_savs(X,pbeta,grp_idx)

% Define Stuff
G = size(unique(grp_idx),1);
M = size(X,2);
T = size(X,1);
grp_size = histc(grp_idx, unique(grp_idx));
pSAVS = NaN(size(pbeta));
n_samples = size(pbeta,2);

% Sparsification Loop
for i = 1:n_samples
for j = 1:G
    ind = find(grp_idx == j);
    l2b = sqrt((pbeta(ind,i)'*pbeta(ind,i)));
    l2x2 = norm(X(ind))^2;
    mu_j = 1/(l2b^2); % was at 1/(l2b^2); was at 1/(l2b^2.5) for the simulations

    if mu_j/(l2x2*l2x2) >= 1
        pSAVS(ind,i) = 0;
    else
     pSAVS(ind,i) = (1-mu_j/mu_j/(l2x2*l2x2))*pbeta(ind,i);
    end
end
end

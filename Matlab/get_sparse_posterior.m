%% Get sparsified posteriors
function [betas_final,pincl_temp] = get_sparse_posterior(data)

post_process = data.post_process;
out = data.out;
grp_idx_temp = data.grp_idx_temp;
midas_prior = data.midas_prior;
sum_grp = data.sum_grp;
Qj = {data.Qj};
Lam_inv_sqr = {data.Lam_inv_sqr};
tin = data.tin;
v = data.v;
MCMC = data.MCMC;
xind = data.xind;
groupall = data.groupall;
pincl_temp = data.pincl_temp;

% Perform group sparsification
if strcmp(post_process,"yes")
[beta_out] = group_savs_orth(out.theta,grp_idx_temp');
else
    beta_out= out.theta;
    betas_final = out.theta;
end

%  Transform back to non-orthogonalised
if ~strcmp(midas_prior,"horseshoe")
betas_final = out.theta;
for j = 1:sum_grp
xind1 = find(grp_idx_temp == j);
betas_final(xind1,:) = Qj{j}*Lam_inv_sqr{j}*beta_out(xind1,:)/sqrt(tin);
end
end

% Variable Selection Info
idx_first_memb = [];
for iii = 1:size(unique(grp_idx_temp),2)
    idx_first_memb = [idx_first_memb;min(find(iii== grp_idx_temp))];
end
pincl_temp(v,unique(groupall(xind))) = (sum(betas_final(idx_first_memb,:)'~=0)/MCMC)' ;

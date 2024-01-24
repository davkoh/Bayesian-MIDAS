function [lamt,nut] = sample_t_errors(ystar,lamt,nut,nub,diag_obs)

T = length(ystar);

        %% Sample t-errors
lamt = 1./gamrnd((nut+1)/2,2./(nut+(1./diag_obs.*speye(T))*(ystar.^2)));

% sample nu
[nut,flag] = sample_nu_MH(lamt,nut,nub);
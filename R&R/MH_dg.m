function [a] = MH_dg(proposal,a,Sigmatau,prior_a,kappa2)
%% Metropolis Hastings algorithm for sampling the a_zetta variable in the double gamma prior
%% Log difference in priors
logdiffprior = 1-prior_a-proposal/prior_a - 1-prior_a-a/prior_a + proposal - a;

%% Compute log marginalised likelihood for proposal
% A
A = (proposal + 0.5)*log(sqrt(proposal*kappa2));
% B
B = (proposal+0.5)*(sqrt(pi)*2) + gammaln(proposal);
% C
C = (a-1/2)*log(abs(Sigmatau));
% D
D = log(besselk(a-1/2,sqrt(a*kappa2)*abs(Sigmatau)));

lognew = A-B+C-D;

%% Compute log marginalised likelihood for previous value
A = (a + 0.5)*log(sqrt(a*kappa2));
% B
B = (a+0.5)*(sqrt(pi)*2) + gammaln(a);
% C
C = (a-1/2)*log(abs(Sigmatau));
% D
D = log(besselk(a-1/2,sqrt(a*kappa2)*abs(Sigmatau)));

logold = A-B+C-D;

logdifflik = lognew - logold;

%% Metropolis Step

if (log(rand)) < logdiffprior + logdifflik
    a = proposal;
else
    a = a;
end

function [ag] = sample_ag_slice(ag,gam,ag_lb,ag_ub,c,d)

log_likelihood = @(x) -log(gamma(x)) + (x-1)*log(gam) -gam + c*log(d) - log(gamma(c)) + (c-1)*log(x) - d*x;;
%log_likelihood = @(x)  c*log(d) - log(gamma(c)) + (c-1)*log(x) - d*x;

ag = uni_slice(ag,log_likelihood,1,inf,ag_lb,ag_ub,[]);

end
function [V] = sample_V2_slice(V,ystar,priorMean,nu_lb,nu_ub,lambda)

log_likelihood = @(x) log(1/sqrt(2*pi*x)) - 1/2 * ((ystar-priorMean)/x)^2 + log(lambda/(2*sqrt(x))) - lambda*sqrt(x);  
%log_likelihood = @(x)  log(lambda/(2*sqrt(x))) - lambda*sqrt(x);  

V = uni_slice(V,log_likelihood,1,inf,nu_lb,nu_ub,[]);

end

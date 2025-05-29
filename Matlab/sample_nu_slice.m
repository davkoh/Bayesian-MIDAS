function [nu] = sample_nu_slice(nu,lam,nu_lb,nu_ub,alpha,beta)

T = size(lam,1);
sum1 = sum(log(lam));
sum2 = sum(1./lam);

log_likelihood = @(x) T*(x/2*log(x/2)-log(gamma(x/2))) - (x/2+1)*sum1 - x/2*sum2 + (alpha-1)*log(x) -beta*x - (log(gamma(alpha)) + alpha*log(beta));

nu = uni_slice(nu,log_likelihood,1,inf,nu_lb,nu_ub,[]);

end



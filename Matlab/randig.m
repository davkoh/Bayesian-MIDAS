function out = randig(theta,chi)
% Generate a draw from the Inverse Gaussian distribution
% Written by Alex Bar-Guy for the RANDRAW toolbox Version 2.0 
%
% The Inverse Gaussian distribution is left skewed distribution whose
% location is set by the mean with the profile determined by the
% scale factor.  The random variable can take a value between zero and
% infinity.  The skewness increases rapidly with decreasing values of
% the scale parameter.
%
%
% pdf(y) = sqrt(chi/(2*pi*y^3)) * exp(-chi./(2*y).*(y/theta-1).^2);
% cdf(y) = normcdf(sqrt(chi./y).*(y/theta-1)) + ...
%            exp(2*chi/theta)*normcdf(sqrt(chi./y).*(-y/theta-1));
%
%   where  normcdf(x) = 0.5*(1+erf(y/sqrt(2))); is the standard normal CDF
%         
% Mean     = theta;
% Variance = theta^3/chi;
% Skewness = sqrt(9*theta/chi);
% Kurtosis = 15*mean/scale;
% Mode = theta/(2*chi)*(sqrt(9*theta^2+4*chi^2)-3*theta);
%
% PARAMETERS:
%  theta - location; (theta>0)
%  chi - scale; (chi>0)
%
% SUPPORT:
%  y,  y>0
%
% CLASS:
%   Continuous skewed distribution
%
% NOTES:
%   1. There are several alternate forms for the PDF, 
%      some of which have more than two parameters
%   2. The Inverse Gaussian distribution is often called the Inverse Normal
%   3. Wald distribution is a special case of The Inverse Gaussian distribution
%      where the mean is a constant with the value one.
%   4. The Inverse Gaussian distribution is a special case of The Generalized
%                    
% Method:
%
% There is an efficient procedure that utilizes a transformation
% yielding two roots.
% If Y is Inverse Gauss random variable, then following to [1]
% we can write:
% V = chi*(Y-theta)^2/(Y*theta^2) ~ Chi-Square(1),
%
% i.e. V is distributed as a chi-square random variable with
% one degree of freedom.
% So it can be simply generated by taking a square of a
% standard normal random number.
% Solving this equation for Y yields two roots:
%
% y1 = theta + 0.5*theta/chi * ( theta*V - sqrt(4*theta*chi*V + ...
%      theta^2*V.^2) );
% and
% y2 = theta^2/y1;
%
% In [2] showed that  Y can be simulated by choosing y1 with probability
% theta/(theta+y1) and y2 with probability 1-theta/(theta+y1)
%
% References:
% [1] Shuster, J. (1968). On the Inverse Gaussian Distribution Function,
%         Journal of the American Statistical Association 63: 1514-1516.
%
% [2] Michael, J.R., Schucany, W.R. and Haas, R.W. (1976).
%     Generating Random Variates Using Transformations with Multiple Roots,
%     The American Statistician 30: 88-90.
%
%

theta=theta';
chi=chi';
sampleSize = length(theta);

chisq1 = randn(sampleSize,1).^2;
out = theta + (0.5*theta./chi).*( theta.*chisq1 - sqrt(4*theta.*chi.*chisq1 + (theta.^2).*(chisq1.^2)) );

l = rand(sampleSize,1) >= theta./(theta+out);
out( l ) = (theta( l ).^2)./out( l );

out=out';

end
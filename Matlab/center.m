function [y,mu,sig] = center(x)

% Function to make your data have mean 0.
% Data in x are Txp, i.e. T time series observations times p variables

mu = nanmean(x);
sig = nanstd(x);

y = (x - repmat(mu,size(x,1),1));

end
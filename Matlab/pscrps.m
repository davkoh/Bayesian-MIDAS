function crps = pscrps(forecasts, obs)
%PSCRPS Calculate Continuous-Ranked Probability Score.
%   CRPS = PSCRPS(FORECAST, OBS) obtains the CRPS for a given vector of
%   observations in OBS, with associated empirical probability forecasts
%   given in the matrix FORECASTS. Each row of FORECASTS contains an
%   empirical distribution for each element of OBS. Output parameter is the
%   CRPS value.
%  
%(cc) Max Little, 2008. This software is licensed under the
%Attribution-Share Alike 2.5 Generic Creative Commons license:
%http://creativecommons.org/licenses/by-sa/2.5/
%If you use this work, please cite:
%Little MA et al. (2008), "Parsimonious Modeling of UK Daily Rainfall for
%Density Forecasting", in Geophysical Research Abstracts, EGU General
%Assembly, Vienna 2008, Volume 10.

[M, S] = size(forecasts);
N = length(obs);

if (M ~= N)
    error('Length of OBS must match number of rows in FORECASTS');
end

crps_individual = zeros(1, N);
for i = 1:N
    crps1 = mean(abs(forecasts(i, :) - obs(i)));
    crps2 = mean(abs(diff(forecasts(i, randperm(S)))));
    crps_individual(i) = crps1 - 0.5*crps2;
end
crps = mean(crps_individual);

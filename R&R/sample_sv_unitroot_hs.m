function [h,Sigmah,h0,lambdah,nuh] = sample_sv_unitroot_hs(Ystar,h,Sigmah,h0,lambdah,nuh,prior_sv)
%% Function Description
% Creates a sample for all parameters related to the SV equation in the
% observation equation using the horseshoe prior on the diagonal of the
% state covariance

% Unpack Priors
Vh = prior_sv.Vh;
T = length(Ystar);

% Sample h
h = SVRW2(Ystar,h,Sigmah,h0);

% sample h0
Kh0 = 1./Sigmah(1)+ 1./Vh;
h0hat = Kh0\(h(1)./Sigmah(1));
h0 = h0hat + chol(Kh0,'lower')'\randn;

% sample Sigmah:
e = (h - [h0;h(1:T-1,:)]);
eta_lamh=1./gamrnd(1,1./(1+1./lambdah));
lambdah=1./gamrnd(1, 1./(1./eta_lamh + (e.^2)./(2.*nuh)));
eta_nuh=1./gamrnd(1,1./(1+1./nuh));
nuh=1./gamrnd((T+1)/2, 1./(1./eta_nuh +0.5*sum(sum(e.^2./lambdah))  ));
Sigmah=(lambdah.*nuh)+1e-10;
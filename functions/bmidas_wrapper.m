%% Wrapper for the BMIDAS Prior Type

function [out] = bmidas_wrapper(input,midas_prior)
% MIDAS prior choice: 
    % gigg
    % horseshoe
    % MAL (Mogliani & Simoni, 2021) model


  if strcmp(midas_prior,"gigg") == 1
      [out] = bmidas_gigg(input);
  elseif strcmp(midas_prior,"horseshoe") == 1
      [out] = bmidas_horseshoe(input);
  elseif strcmp(midas_prior,"MAL") ==1
   % Display the warning message
    warning('For MAL prior, the BMIDAS model does not contain a trend or stochastic volatility.');
      It.nsave = input.samples;
      It.nburn = input.burnin;
      It.nthin = 1;
      Xv = input.X;
      y  = input.Y;
      Xf = Xv(end,:);
      grp_idx = input.grp_idx';
      [out] = BMIDAS_AGLasso_SS_independent(Xv,y,Xf,It,grp_idx);
      % Match downstream conventions: theta is K x MCMC.
      out.theta = out.beta';
      out.alpha = zeros(input.samples,1); % intercept already included in nowcasts
  else
        error('Unsupported midas_prior value: %s. Use "gigg", "horseshoe", or "MAL".', midas_prior);
    end
  end

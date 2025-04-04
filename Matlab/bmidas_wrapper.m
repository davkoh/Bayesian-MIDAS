%% Wrapper for the BMIDAS Prior Type

function [out] = bmidas_wrapper(input,midas_prior)
% MIDAS prior choice: 
    % gigg
    % horseshoe
    % MAL (Mogliani & Simoni, 2022) model


  if strcmp(midas_prior,"gigg") == 1
      [out] = bmidas_gigg(input);
  elseif strcmp(midas_prior,"horseshoe") == 1
      [out] = bmidas_horseshoe(input);
  elseif strcmp(midas_prior,"MAL") ==1
   % Display the warning message
    warning('For MAL prior, the BMIDAS model does not contain a trend or stochastic volatility.');
      [out] = bmidas_mal(input);
  else
        error('Unsupported midas_prior value: %s. Use "gigg", "horseshoe", or "MAL".', midas_prior);
    end
  end

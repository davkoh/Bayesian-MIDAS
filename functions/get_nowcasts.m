function  [crpsv,predv,trendv,svv,sv_trendv,cycpredv] = get_nowcasts(data)

betas_final = data.betas_final;
out = data.out;
sv_obs_type =data.sv_obs_type;
midas_type = data.midas_type;
Xm = data.Xm;
grp_idx = data.grp_idx;
poly = data.poly;
tperiod = data.tperiod;
xind = data.xind;
v = data.v;
tail_type = data.tail_type;
trend_type = data.trend_type;
tin = data.tin;
MCMC = data.MCMC;
yf = data.yf;
crpsv = data.crpsv;
predv= data.predv;
midas_prior = data.midas_prior;
cycpredv = data.cycpredv;
trendv = data.trendv;
svv = data.svv;
sv_trendv = data.trendv;


if strcmp(midas_prior,"gigg") || strcmp(midas_prior,"horseshoe") 

% Retrieve parameters
betas_final = betas_final';
hout = out.h;
thetaout = out.state_params;
nuout = out.nu;
if ~strcmp(trend_type,"none")
    gout = out.g;
    tauout = out.tau;
end
if strcmp(sv_obs_type,"none")
    sigma2 = out.sigma2;
end

ypredtt = []; % local storage
crps_temp = []; % local storage for crps values
cycpredtt = []; % local storage
trendtt = []; % local storage for trend 
svtt = []; % local storage for trend 
sv_trendtt = []; % local storage for trend 

if midas_type == "almon"
% Get out of sample Almon data
if v < 2
[Xv,~] = midas_dat_r2_final(Xm(1:tin-3+tperiod,xind),grp_idx,poly);
else
    [Xv,~] = midas_dat_r2_final(Xm(1:tin+tperiod,xind),grp_idx,poly);
end
Xv = (Xv);
Xv = Xv(end,:);

else
    if v < 2
    Xv = Xm(tin-3+tperiod,xind);
    else
        Xv = Xm(tin+tperiod,xind);
    end
end

t_cont = 1;
sv_cont = 1;
trend_cont = 0;

% Monte Carlo Integration for predictive distribution
for j = 1:(MCMC)
  
    % Sample g_{t+1} and tau_{t+1} (only when trend is active)
    if ~strcmp(trend_type,"none")
        if v<2
          g_temp = gout(end,j) + randn*thetaout(j,2);
          tau_temp = tauout(end,j) + exp(0.5*g_temp)*randn;
          g_temp = g_temp + randn*thetaout(j,2);
        else
            g_temp = gout(end,j) + randn*thetaout(j,2);
        end

        % Sample tau_{t+1}
        if v<2
          tau_temp = tau_temp  + exp(0.5*g_temp)*randn;
        else
            tau_temp = tauout(end,j) + exp(0.5*g_temp)*randn;
        end
    end

     % Sample h_{t+1} (only when observation SV is active)
     if ~strcmp(sv_obs_type,"none")
         if v<2
          h_temp = hout(end,j) + randn*thetaout(j,1);
          h_temp = h_temp + randn*thetaout(j,1);
         else
             h_temp = hout(end,j) + randn*thetaout(j,1);
         end
     end

     % Sample y_{t+1}

     if strcmp(tail_type,"terr")
         t_cont = trnd(nuout(j));
     else
         t_cont = randn;
     end

     if ~strcmp(sv_obs_type,"none")
         sv_cont = exp(0.5*h_temp); % add a line for the scale of linear regression
     else
         sv_cont = sqrt(sigma2(j));
     end
     if ~strcmp(trend_type,"none")
         trend_cont = tau_temp;
     else 
         trend_cont = out.alpha(j);
     end

      y_temp = trend_cont + Xv*betas_final(j,:)' + sv_cont*t_cont;

      ypredtt= [ypredtt y_temp];
      cycpredtt= [cycpredtt Xv*betas_final(j,:)'];
      if ~strcmp(trend_type,"none")
          trendtt = [trendtt tau_temp];
          sv_trendtt = [sv_trendtt g_temp];
      end
      if ~strcmp(sv_obs_type,"none")
          svtt = [svtt h_temp];
      end
      if strcmp(tail_type,"terr")
      crps_temp = [crps_temp; crps_t(yf(1,1),nuout(j),trend_cont + Xv*betas_final(j,:)',sv_cont)];
      else
          crps_temp = [crps_temp; crps_t(yf(1,1),200,trend_cont + Xv*betas_final(j,:)',sv_cont)];
      end
end

end

if strcmp(midas_prior, "MAL")
betas_final = betas_final;
sig2 = out.sigma2;
intercept = out.beta0;

ypredtt = []; % local storage
crps_temp = []; % local storage for crps values

if midas_type == "almon"
% Get out of sample Almon data
if v < 2
[Xv,~] = midas_dat_r2_final(Xm(1:tin-3+tperiod,xind),grp_idx,poly);
else
    [Xv,~] = midas_dat_r2_final(Xm(1:tin+tperiod,xind),grp_idx,poly);
end
Xv = (Xv);
Xv = Xv(end,:);
else
        if v < 2
    Xv = Xm(tin-3+tperiod,xind);
    else
        Xv = Xm(tin+tperiod,xind);
    end
end

t_cont = 1;
sv_cont = 1;
trend_cont = 0;

% Monte Carlo Integration for predictive distribution
for j = 1:(MCMC)
  

  
         t_cont = randn;
  

         sv_cont = sqrt(sig2(j));
     
         trend_cont = intercept(j);
     

      y_temp = trend_cont + Xv*betas_final(j,:)' + sv_cont*t_cont;

      ypredtt= [ypredtt y_temp];
      
          crps_temp = [crps_temp; crps_t(yf(1,1),200,trend_cont + Xv*betas_final(j,:)',sv_cont)];
      
 end

end

%%%%%%%%%%%%  Save prediction results by nowcast period %%%%%%%%%%%%%%%
crpsv = [crpsv;mean(crps_temp)];
predv = [predv;ypredtt];

if ~strcmp(midas_prior,"MAL")
    cycpredv = [cycpredv;cycpredtt]; 
    if ~strcmp(trend_type,"none")
        trendv = [trendv;trendtt];
        sv_trendv = [sv_trendv;sv_trendtt];
    end
    if ~strcmp(sv_obs_type,"none")
        svv = [svv;svtt];
    end
end

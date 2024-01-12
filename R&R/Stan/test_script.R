## Test Script using Stan for the HS - LT - SV model

library(cmdstanr)
library(resample)
library(loo)
library(posterior)
library(scoringRules)
library(doParallel)

rm(list = ls())

setwd("D:/Github/Bayesian-MIDAS/R&R")

# Generate Data
T = 200;
sd_tau = 1;
t_0 = 0;
sd_y = 1;
h_0 = 0;
sd_h = 1;


tau = array(0,c(T,1));
y = array(0,c(T,1));
h = array(0,c(T,1));


for (t in 1:T){
if (t == 1){
tau[1] <- t_0 + sd_tau*rnorm(1)
h[1] <- h_0 + sd_h*rnorm(1)
} else{
tau[t] <- tau[t-1] + sd_tau*rnorm(1);
h[t] <- 0.95*h[t-1] + sd_h*rnorm(1);
}


if (t == 100){
tau[t] <- tau[t] + 10;
h[t] <- h[t];
y[t] <- tau[t] + exp(1/2*h[t])*rnorm(1);
} else {
  y[t] <- tau[t] + exp(1/2*h[t])*rnorm(1);
}
  }


tau_true = tau;
h_true = h;

ts.plot(y)
ts.plot(h_true)
ts.plot(tau_true)

# Compile and Estimate Model

dat <- list(
  T = T,
  y = as.vector(y)
)

mod <- cmdstan_model("Stan/lt_sv_hs_trans.stan")

fit <- mod$sample(data = dat,
                  seed = 123,
                  chains = 4,
                  parallel_chains = 4,
                  iter_warmup = 500,
                  iter_sampling = 500,
                  adapt_delta = 0.95,
                  max_treedepth = 15
                  )

fit$summary(variables = c("nu_tau"))

x <- fit$draws(variables = c("phi"), format = "matrix")

bayesplot::mcmc_trace(x, pars = c("phi"))

ptau <- colMeans(fit$draws(variables = c("tau"),format = "matrix"))

ts.plot(ptau)
ts.plot(tau_true)


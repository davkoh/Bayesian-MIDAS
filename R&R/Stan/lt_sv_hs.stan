data{
  int<lower=1> T;
  vector[T] y;
}

parameters{
  real h0;
  real tau0;
  vector[T] z_tau;
  vector[T] z_h;
  real<lower=0> nu_h;
  real<lower=0> nu_tau;
//  vector<lower=0>[T] lam_h;
//  vector<lower=0>[T] lam_tau;
}

transformed parameters{
  real<lower=0> sig_tau = nu_tau;//nu_tau * lam_tau;
  real<lower=0> sig_h = nu_h;//nu_h * lam_h;
  vector[T] tau;
  vector[T] h;
  
  h[1] = h0 + sig_h*z_h[1];
  
  for (t in 2:T){
    h[t] = h[t-1] + sig_h*z_h[t];
  }

  tau[1] = tau0 + sig_tau*z_tau[1];
  
    for (t in 2:T){
    tau[t] = tau[t-1] + sig_tau*z_tau[t];
  }
  
}

model{
  // Priors
  nu_h ~ normal(0,1);
  nu_tau ~ normal(0,1);
//  lam_h ~ normal(0,1);
//  lam_tau ~ normal(0,1);
  // Likelihood
  h0 ~ normal(0,10);
  tau0 ~ normal(0,10);
  y ~ normal(tau,exp(h / 2));
}


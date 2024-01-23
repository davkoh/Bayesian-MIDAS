data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real mu_h;                     // mean log volatility
  real mu_g;
  real mu_t;
  real<lower=-1,upper=1> phi_h;  // persistence of volatility
  real<lower=-1,upper=1> phi_g;  // persistence of volatility
  vector[T] z_h;
  vector[T] z_g;
  vector[T] z_t;
  real<lower=0> nu_g;
  real<lower=0> nu_h;
  vector<lower=0>[T] lam_h;
  
}

transformed parameters{
vector<lower=0>[T] sig_h = nu_h * lam_h;
real<lower=0> sig_g = nu_g;
vector[T] h = z_h .* sig_h;  // now h ~ normal(0, sigma)
vector[T] g = z_g * sig_g;
vector[T] tau = z_t .* exp(g/2);

// Assemble h
  h[1] /= sqrt(1 - phi_h * phi_h);  // rescale h[1]
  h += mu_h;
  for (t in 2:T) {
    h[t] += phi_h * (h[t-1] - mu_h);
}

// Assemble g

  g[1] /= sqrt(1 - phi_g * phi_g);  // rescale g[1]
  g += mu_g;
  for (t in 2:T) {
    g[t] += phi_g * (g[t-1] - mu_g);
}

// Assemble tau

  tau += mu_t;
  for (t in 2:T) {
    tau[t] += (tau[t-1] - mu_t) ;
}

}

model {
  // AR Paramters 
  phi_h ~ normal(0,2);
  phi_g ~ normal(0, 2);
  // Constants
  mu_h ~ normal(0,2);
  mu_t ~ normal(0,2);
  mu_g ~ normal(0,2);
  // Hyperparameters
  nu_g ~ student_t(4,0,2);
  nu_h ~ student_t(4,0,2);
  lam_h ~ student_t(4,0,2);
  // Non-centred terms
  z_h ~ std_normal();
  z_t ~ std_normal();
  z_g ~ std_normal();
  // Likelihood
  y ~ normal(tau, exp(h / 2));
}

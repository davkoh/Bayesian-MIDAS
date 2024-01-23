data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real h0;
  real g0;
  real t0;
  vector[T] z_h;
  vector[T] z_g;
  vector[T] z_t;
  real<lower=0> nu_g;
  real<lower=0> nu_h;
  vector<lower=0>[T] lam_h;
  vector<lower=0>[T] lam_g; 
}

transformed parameters{
vector<lower=0>[T] sig_h = nu_h * lam_h;
vector<lower=0>[T] sig_g = nu_g * lam_g;
vector[T] h;  // now h ~ normal(0, sigma)
vector[T] g;
vector[T] tau;

// Assemble h
  h[1] = h0 + z_h[1]* sig_h[1];  // rescale h[1]
  for (t in 2:T) {
    h[t] = h[t-1] + z_h[t]*sig_h[t];
}

// Assemble g
  g[1] = g0 + z_g[1]* sig_g[1];  // rescale h[1]
  for (t in 2:T) {
    g[t] = g[t-1] + z_g[t]*sig_g[t];
}

// Assemble tau

  tau[1] = t0 + z_t[1]*exp(g[1]/2);
  for (t in 2:T) {
    tau[t] = tau[t-1] + z_t[t]*exp(g[t]/2);
}

}

model {
  // Initial conditions
  t0 ~ normal(0,10);
  g0 ~ normal(0,10);
  h0 ~ normal(0,10);
  // Hyperparameters
  lam_g ~ student_t(4,0,2);
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

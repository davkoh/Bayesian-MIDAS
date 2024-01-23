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
  real<lower=0> sig_h;
  real<lower=0> sig_g;
}

transformed parameters{
vector[T] h;  // now h ~ normal(0, sigma)
vector[T] g;
vector[T] tau;

// Assemble h
  h[1] = h0 + z_h[1]* sig_h;  // rescale h[1]
  for (t in 2:T) {
    h[t] = h[t-1] + z_h[t]*sig_h;
}

// Assemble g
  g[1] = g0 + z_g[1]* sig_g;  // rescale h[1]
  for (t in 2:T) {
    g[t] = g[t-1] + z_g[t]*sig_g;
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
  sig_g ~ normal(0,1);
  sig_h ~ normal(0,1);
  // Non-centred terms
  z_h ~ std_normal();
  z_t ~ std_normal();
  z_g ~ std_normal();
  // Likelihood
  y ~ normal(tau, exp(h / 2));
}

data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real mu;                     // mean log volatility
  real mu_tau;
  real<lower=-1,upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  vector[T] z_h;
  vector[T] z_tau;
  //vector<lower=0>[T] sig_tau;
  real<lower=0> nu_tau;
  vector<lower=0>[T] lam_tau; 
}

transformed parameters{
vector<lower=0>[T] sig_tau = nu_tau * lam_tau;
vector[T] h = z_h * sigma;  // now h ~ normal(0, sigma)
vector[T] tau = z_tau .* sig_tau;
  h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
  h += mu;
  for (t in 2:T) {
    h[t] += phi * (h[t-1] - mu);
}


  tau += mu_tau;
  for (t in 2:T) {
    tau[t] += (tau[t-1] - mu_tau);
}

}

model {
  phi ~ normal(0, 2);
  sigma ~ student_t(4,0, 2);
  lam_tau ~ student_t(4,0,2);
  nu_tau ~ student_t(4,0,2);
  mu ~ normal(0, 10);
  mu_tau ~ normal(0,10);
  z_h ~ std_normal();
  z_tau ~ std_normal();
  y ~ normal(tau, exp(h / 2));
}


data{
  int<lower=0> T;
  vector[T] y;
}

parameters{
  real<lower=0> sigma_y;
  real<lower=-1,upper=1> phi;
  real<lower=0> nu; // global shrinkage
  vector<lower=0,upper=1>[T] tilde_R; // used for transformation to beta prime
  //real<lower=0,upper=1> tilde_R2;
  vector[T] z_tau;
  
}

transformed parameters{
  vector[T] zi =  log( tilde_R ./ (1-tilde_R));
  real mu = log(nu);
  vector[T] h;
  vector[T] tau;
  
  h[1] = mu + zi[1];
  for (t in 2:T){
        h[t] = mu + phi * (h[t - 1] -  mu) +  zi[t];
  }
  
  tau[1] = z_tau[1]*exp(1.0/2*h[1]);
    for (t in 2:T){
        tau[t] = tau[t-1] + exp(1.0/2*h[t])*z_tau[t];
  }
  
}

model{
  y ~ normal(tau,sigma_y);
  nu ~ student_t(1,0,sigma_y/sqrt(T));
  tilde_R ~ beta(0.1,0.5);
  //tilde_R2 ~ beta(alpha,beta);
  (phi +1)/2 ~ beta(10,2);
  sigma_y ~ student_t(3,0,sd(y));
  z_tau ~ std_normal();
}

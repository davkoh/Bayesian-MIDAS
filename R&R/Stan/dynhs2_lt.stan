
data{
  int<lower=0> T;
  vector[T] y;
  real<lower=0> alpha; // first parameter of the beta
  real<lower=0> beta; // second parameter of the beta
}

parameters{
  real<lower=0> sigma_y;
  real<lower=-1,upper=1> phi;
  real<lower=0> nu; // global shrinkage
  //vector<lower=0,upper=1>[T] tilde_R; // used for transformation to beta prime
  //real<lower=0,upper=1> tilde_R2;
  vector[T] z_tau;
  vector<lower=0>[T] lam;
  
  
}

transformed parameters{
  vector[T] zi =  log( lam);
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
  lam ~ student_t(4,0,2);
  //tilde_R2 ~ beta(alpha,beta);
  (phi +1)/2 ~ beta(10,2);
  sigma_y ~ student_t(3,0,sd(y));
  z_tau ~ std_normal();
}

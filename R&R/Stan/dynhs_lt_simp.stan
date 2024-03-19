
data{
  int<lower=0> T;
  vector[T] y;
  real<lower=0> alpha; // first parameter of the beta
  real<lower=0> beta; // second parameter of the beta
}

parameters{
  real<lower=0> sigma_y;
  //real<lower=-1,upper=1> phi;
  real<lower=0> nu; // global shrinkage
  vector<lower=0>[T] lam; // used for transformation to beta prime
  //real<lower=0,upper=1> tilde_R2;
  vector[T] z_tau;
  
  
}

transformed parameters{
  vector[T] h = lam * nu;
  vector[T] tau;
  

  tau[1] = z_tau[1]*h[1];
    for (t in 2:T){
        tau[t] = tau[t-1] + h[t]*z_tau[t];
  }
  
}

model{
  y ~ normal(tau,sigma_y);
  nu ~ student_t(1,0,1);
  lam ~ student_t(1,0,1);
  sigma_y ~ student_t(3,0,sd(y));
  z_tau ~ std_normal();
}

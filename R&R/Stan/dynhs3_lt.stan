
data{
  int<lower=0> T;
  vector[T] y;
  real<lower=0> alpha; // first parameter of the beta
  real<lower=0> beta; // second parameter of the beta
}

parameters{
  real<lower=-1,upper=1> phi;
  real<lower=0> nu; // global shrinkage
  //vector<lower=0,upper=1>[T] tilde_R; // used for transformation to beta prime
  //real<lower=0,upper=1> tilde_R2;
  vector[T] z_tau;
  vector<lower=0>[T] lam;
  real<lower=0> sigma_vol;
  real<lower=-1,upper=1> phi_vol;
  vector[T] vol_std;
  real mu_vol;
  
  
}

transformed parameters{
  vector[T] zi =  log( lam);
  real mu = log(nu);
  vector[T] vol = vol_std * sigma_vol;  // now h ~ normal(0, sigma)
  vector[T] h;
  vector[T] tau;
  
   vol[1] /= sqrt(1 - phi_vol * phi_vol);  // rescale h[1]
  vol += mu_vol;
  for (t in 2:T) {
    vol[t] += phi_vol * (vol[t - 1] - mu_vol);
  }
  
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
  y ~ normal(tau,exp(vol/2));
  nu ~ student_t(1,0,mean(exp(vol/2))/sqrt(T));
  lam ~ student_t(4,0,2);
  //tilde_R2 ~ beta(alpha,beta);
  (phi +1)/2 ~ beta(10,2);
  z_tau ~ std_normal();
  
  // Vol equation
  vol_std ~ std_normal();
  (phi_vol +1)/2 ~ beta(10,2);
  sigma_vol ~ normal(0,3);
  mu_vol ~ normal(0, 10);

}

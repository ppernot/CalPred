data {
  int N;
  vector[N] x;
  vector[N] y;
  vector[N] uy;
  int N2;
  vector[N2] x_new;
  int<lower=0, upper=1> inadequacy;
}
parameters {
  real eps;
  real sig;
  real<lower=-4> log_sigma;
}
transformed parameters{
  real sigma;
  vector[N] vy;
  sigma = 10^log_sigma;
  vy = phys_mod_vec(x,eps,sig,inadequacy);
}
model {
  y ~ normal(vy,sigma); 
}
generated quantities{
  vector[N]  resid;
  vector[N]  y_pred_cont;
  vector[N2] y_conf;
  vector[N2] y_pred;
  real       br;
  real       up;
  
  # Residuals
  resid = y - vy;

  # Birge ratio
  {
    vector[N]  Vm1;
    for (i in 1:N)
      Vm1[i] = 1/sigma^2;
    br = quad_form(diag_matrix(Vm1),resid) / (N-3);
  }
  
  # Prediction uncert stat 
  for (i in 1:N)
    y_pred_cont[i] = normal_rng(vy[i],sigma);

  # Predicted data
  y_conf = phys_mod_vec(x_new,eps,sig,inadequacy);
  for (i in 1:N2)
    y_pred[i] = normal_rng(y_conf[i],sigma);
}

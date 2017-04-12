data {
  int N;
  vector[N] x;
  vector[N] y;
  vector[N] uy;
  int N2;
  vector[N2] x_new;
  int<lower=0, upper=1> inadequacy;
  int<lower=0, upper=1> shift; # Affects prediction
}
parameters {
  real eps;
  real sig;
  real<lower=-4> log_sigma;
}
transformed parameters{
  vector[N] mu_M;
  vector[N] V_D;
  real sigma;
  sigma = 10^log_sigma;
  
  mu_M = phys_mod_vec(x,eps,sig,inadequacy);
  for (i in 1:N)
      V_D[i] = uy[i]^2 + sigma^2;

}
model {
  for (i in 1:N)
    y[i] ~ normal(mu_M[i], V_D[i]^0.5);
}
generated quantities{
  vector[N]  resid;
  vector[N]  y_pred_cont;
  vector[N2] y_conf;
  vector[N2] y_pred;
  real       br;
  real       vy;
  real       uyp;

  # Residuals
  resid = y - mu_M;

  # Birge ratio
  br = quad_form(inverse(diag_matrix(V_D)),resid) / (N-3);
  
  # Pres stat
  for (i in 1:N)
    y_pred_cont[i] = normal_rng(mu_M[i],V_D[i]^0.5);

  # Predicted data
  for (i in 1:N2) {
    vy  = phys_mod(x_new[i],eps,sig,inadequacy);
    uyp = fmax(min(uy),vy * mean(uy ./ y));
    if(shift ==0) {
      y_conf[i] = normal_rng(vy,sigma);
      y_pred[i] = normal_rng(y_conf[i],uyp);
    } else {
      y_conf[i] = vy;
      y_pred[i] = normal_rng(y_conf[i], (uyp^2+sigma^2)^0.5);
    }
  }
}

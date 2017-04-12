data {
  int N;
  vector[N] x;
  vector[N] y;
  vector[N] uy;
  int N2;
  vector[N2] x_new;
  int<lower=0, upper=1> inadequacy;
  real sigma;
}
parameters {
  real eps;
  real sig;
}
model {
  y ~ normal(phys_mod_vec(x,eps,sig,inadequacy),sigma); 
}
generated quantities{
  vector[N]  resid;
  vector[N2] y_conf;
  vector[N2] y_pred;
  real       br;

  # Residuals
  resid = y - phys_mod_vec(x,eps,sig,inadequacy);

  # Birge ratio
  {
    vector[N]  Vm1;
    for (i in 1:N)
      Vm1[i] = 1/sigma^2;
    br = quad_form(diag_matrix(Vm1),resid) / (N-2);
  }
  
  # Predicted data
  y_conf = phys_mod_vec(x_new,eps,sig,inadequacy);
  y_pred = y_conf;
}

data {
  int N;
  vector[N] x;
  vector[N] y;
  vector[N] uy;
  int tag[N];
  int N2;
  vector[N2] x_new;
  int<lower=0, upper=1> inadequacy;
}
transformed data {
  int<lower=1> n_grp;  // nb groups
  n_grp = max(tag);
}
parameters {
  real eps;
  real sig;
  real<lower= -4, upper=0> log_tau;
}
transformed parameters {
  real tau;
  vector[N] mu_M;
  cov_matrix[N] U;
  
  tau = 10^log_tau;
  
  # Data cov matrix
  for (k in 1:(N-1)) {
    for (l in (k+1):N) {
      U[k,l] = 0;
      if(tag[k]==tag[l]) U[k,l] = tau^2 ;
      U[l,k] = U[k,l];
    }
    U[k,k] = tau^2 + uy[k]^2;
  }
  U[N,N] = tau^2 + uy[N]^2;
  
  mu_M = phys_mod_vec(x,eps,sig,inadequacy);
  
}
model {
  y ~ multi_normal(mu_M, U);
}
generated quantities{
  vector[N] resid;
  vector[N] y_pred_cont;
  vector[N2] y_conf;
  vector[N2] y_pred;
  real vy;
  real uyp;
  real br;

  resid = y - mu_M; 

  br = quad_form(inverse(U),resid) / (N-3);
  
  for (i in 1:N)
    y_pred_cont[i] = normal_rng(mu_M[i],U[i,i]^0.5);

  for (i in 1:N2) {
    vy  = phys_mod(x_new[i],eps,sig,inadequacy);
    uyp = fmax(min(uy),vy * mean(uy ./ y));
    y_conf[i] = vy;
    y_pred[i] = normal_rng(y_conf[i], (uyp^2+tau^2)^0.5);
  }
}

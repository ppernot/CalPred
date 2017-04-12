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
  real<lower=-4> log_u_eps;
  real<lower=-4> log_u_sig;
  real<lower=-4, upper=0> log_tau;
  vector[n_grp] shift;
  cholesky_factor_corr[2] L;
}
transformed parameters {
  cov_matrix[2] V;
  corr_matrix[2] C;
  vector[2] mu;
  vector[2] g;
  vector[N] mu_M;
  vector[N] V_T;
  vector[n_grp] shiftc;
  vector[2] sigma;
  real u_M2;
  real tau;
  
  shiftc = shift - mean(shift);
  
  sigma[1] = 10^log_u_eps; sigma[2] = 10^log_u_sig;
  tau = 10^log_tau;
  
  C = L * L';
  V = diag_matrix(sigma) * C  * diag_matrix(sigma);

  mu[1] = eps; mu[2] = sig;
  
  for (i in 1:N) {
    mu_M[i] = phys_mod(x[i],eps,sig,inadequacy) + shiftc[tag[i]];
    g = grad_phys_mod(x[i],eps,sig,inadequacy);
    u_M2 = g' * V * g;
    V_T[i] = uy[i]^2 + u_M2;
  }

}
model {
  L ~ lkj_corr_cholesky(1);
  
  to_vector(shift) ~ normal(0,tau);
  
  for (i in 1:N) 
    y[i] ~ normal(mu_M[i], V_T[i]^0.5);
  
}
generated quantities{
  vector[N] resid;
  vector[N2] y_conf;
  vector[N2] y_pred;
  vector[N] y_pred_cont;
  real vy;
  real uyp;
  vector[2] mup;
  real br;
  real u_eps;
  real u_sig;
  real rho;

  u_eps = sigma[1]; u_sig = sigma[2]; rho = C[1,2];

  resid = y - mu_M;
  
  br = quad_form(inverse(diag_matrix(V_T)),resid)  / (N-n_grp-5);

  for (i in 1:N)
    y_pred_cont[i] = normal_rng(mu_M[i],V_T[i]^0.5);
    
  mup = multi_normal_rng(mu,V);
  for (i in 1:N2) {
    vy  = phys_mod(x_new[i],mup[1],mup[2],inadequacy);
    uyp = fmax(min(uy),vy * mean(uy ./ y));
    y_conf[i] = vy;
    y_pred[i] = normal_rng(vy,(uyp^2 + tau^2)^0.5);
  }
}

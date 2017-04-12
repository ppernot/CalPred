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
  real<lower= -4> log_u_eps;
  real<lower= -4> log_u_sig;
  real<lower= -4, upper=0> log_tau;
  cholesky_factor_corr[2] L_Omega;
  matrix[2,n_grp] z;
}
transformed parameters {
  real u_eps;
  real u_sig;
  real tau;
  vector[2] tauMu;
  matrix[2,n_grp] mu; 
  matrix[2,n_grp] mut;
  cov_matrix[N] U;
  vector[N] mu_M;
  
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
  
  u_eps = 10^log_u_eps;
  tauMu[1] = u_eps;
  u_sig = 10^log_u_sig;
  tauMu[2] = u_sig;
  
  for (i in 1:n_grp) {
    mut[1,i] = eps;
    mut[2,i] = sig;
  }
  mu = mut + diag_pre_multiply(tauMu,L_Omega) * z; 
 
  for (i in 1:N)
    mu_M[i] = phys_mod(x[i],mu[1,tag[i]],mu[2,tag[i]],inadequacy); 
  
}
model {
  log_tau ~ normal(-4,1);
  
  to_vector(z) ~ normal(0,1);

  L_Omega ~ lkj_corr_cholesky(1);

  y ~  multi_normal(mu_M, U);
   
}
generated quantities{
  vector[N] resid;
  vector[N] y_pred_cont;
  vector[N2] y_conf;
  vector[N2] y_pred;
  vector[2] mup;
  vector[2] zp;
  corr_matrix[2] Cor;
  real vy;
  real uyp;
  real br;
  real rho;
  
  Cor = L_Omega * L_Omega';
  rho = Cor[1,2];
  
  resid = y - mu_M; 

  br = quad_form(inverse(U),resid) / (N-2*n_grp-1);
  
  // for (i in 1:N)
  //   y_pred_cont[i] = normal_rng(mu_M[i], U[i,i]^0.5);

  for (k in 1:2)
    zp[k] = normal_rng(0,1);
  mup = mut[,1] + diag_pre_multiply(tauMu,L_Omega) * zp;  

  for (i in 1:N) {
    vy  = phys_mod(x[i],mup[1],mup[2],inadequacy);
    y_pred_cont[i] = normal_rng(vy, U[i,i]^0.5);
  }
  
  for (i in 1:N2) {
    vy  = phys_mod(x_new[i],mup[1],mup[2],inadequacy);
    uyp = fmax(min(uy),vy * mean(uy ./ y));
    y_conf[i] = vy;
    y_pred[i] = normal_rng(vy, (uyp^2 + tau^2)^0.5);
  }
}

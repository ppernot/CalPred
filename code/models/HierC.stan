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
  real<lower=-4> log_sigma;
  cholesky_factor_corr[2] L_Omega;
  matrix[2,n_grp] z;
}
transformed parameters {
  vector[2] tauMu;
  matrix[2,n_grp] mu; 
  matrix[2,n_grp] mut;
  real u_eps;
  real u_sig;
  real sigma;
  
  u_eps = 10^log_u_eps;
  u_sig = 10^log_u_sig;
  sigma = 10^log_sigma;
  
  tauMu[1] = u_eps;
  tauMu[2] = u_sig;
  
  for (i in 1:n_grp) {
    mut[1,i] = eps;
    mut[2,i] = sig;
  }
  mu = mut + diag_pre_multiply(tauMu,L_Omega) * z; 
}
model {

  to_vector(z) ~ normal(0,1);

  L_Omega ~ lkj_corr_cholesky(1);
    
  for (i in 1:N) 
    y[i] ~ normal(phys_mod(x[i],mu[1,tag[i]],mu[2,tag[i]],inadequacy),
                  (uy[i]^2 + sigma^2)^0.5 );
  
}
generated quantities{
  vector[N] resid;
  vector[N] mu_M;
  vector[N] y_pred_cont;
  vector[N2] y_conf;
  vector[N2] y_pred;
  vector[2] mup;
  vector[2] zp;
  real vy;
  real uyp;
  vector[N] Vm1;
  real br;
  corr_matrix[2] Cor;
  real rho;
  
  Cor = L_Omega * L_Omega';
  rho = Cor[1,2];

  // Local residuals
  for (i in 1:N) 
    resid[i] = y[i] - phys_mod(x[i],mu[1,tag[i]],mu[2,tag[i]],inadequacy);
  for (i in 1:N) 
    Vm1[i] = 1/(uy[i]^2 + sigma^2);

  br = quad_form(diag_matrix(Vm1),resid) / (N-2*n_grp-1);

  for (k in 1:2)
    zp[k] = normal_rng(0,1);
  mup = mut[,1] + (diag_pre_multiply(tauMu,L_Omega) * zp);  

  mu_M  = phys_mod_vec(x,mup[1],mup[2],inadequacy);
  for (i in 1:N)
    y_pred_cont[i] = normal_rng(mu_M[i],Vm1[i]^-0.5);
    
  for (i in 1:N2) {
    if(inadequacy == 1) { # Params compensate for model inadequacy 
      vy  = phys_mod(x_new[i],mup[1],mup[2],inadequacy);
      uyp = fmax(min(uy),vy * mean(uy ./ y));
      y_conf[i] = vy ; 
      y_pred[i] = normal_rng(y_conf[i], (uyp^2 + sigma^2)^0.5);
    } else {              # Params describe data inconsistency
      vy  = phys_mod(x_new[i],eps,sig,inadequacy);
      uyp = fmax(min(uy),vy * mean(uy ./ y));
      y_conf[i] = normal_rng(vy,sigma); 
      y_pred[i] = normal_rng(phys_mod(x_new[i],mup[1],mup[2],inadequacy), uyp);
    }
 
  }
}

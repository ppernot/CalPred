data {
  int N;
  vector[N] x;
  vector[N] y;
  vector[N] uy;
  int N2;
  vector[N2] x_new;
  int<lower=0, upper=1> inadequacy;
  vector[2] primu;
  cov_matrix[2] covprimu;
}
transformed data {
  real v_scale;
  vector[N] x_scaled;
  vector[N2] x_new_scaled;
  
  // Rescale x and variance for GP
  x_scaled = (x-min(x))/(max(x)-min(x));
  x_new_scaled = (x_new-min(x))/(max(x)-min(x));
  v_scale = max(uy)^2;
  
}
parameters {
  vector[2] mu0;
  real<lower=0> eta_sq;
  real<lower=0> inv_rho_sq;
}
transformed parameters {
  real rho_sq;
  vector[N] mu_M;
  cov_matrix[N] V;
 
  rho_sq = inv(inv_rho_sq);
  
  for (i in 1:N) 
    for (j in 1:N) 
      V[i,j] = eta_sq * v_scale * exp(-rho_sq * pow(x_scaled[i] - x_scaled[j],2)) 
             + (i==j ? uy[i]^2 : 0.0);
               
  mu_M = phys_mod_vec(x,mu0[1],mu0[2],inadequacy);

}
model {
  mu0 ~ multi_normal(primu,covprimu);
  
  eta_sq ~ normal(0,1);
  inv_rho_sq ~ normal(0,1);
  
  y ~ multi_normal(mu_M, V);
}
generated quantities{
  vector[N] mu_M0;
  vector[N] resid;
  vector[N] y_pred_cont;
  vector[N2] y_m;
  vector[N2] y_gp;
  vector[N2] y_conf;
  vector[N2] y_pred;
  real vy;
  real uyp;
  real eps;
  real sig;
  matrix[N,N]   Vm1;
  matrix[N,N]   Kr;
  matrix[N,N]   Kr2;
  matrix[N,N2]  K;
  matrix[N2,N]  K2;
  real br;
  real neff;
  
  vector[N] log_lik;
  
  for (i in 1:N)
    log_lik[i] = normal_log(y[i],mu_M[i],V[i,i]^0.5);


  Vm1 = inverse(V);

  for (i in 1:N)
    for (j in 1:N)
      Kr[i,j] = eta_sq * v_scale 
              * exp(-rho_sq *pow(x_scaled[i] - x_scaled[j],2));
  Kr2 = Kr' * Vm1;

  for (i in 1:N)
    for (j in 1:N2)
      K[i,j] = eta_sq * v_scale 
             * exp(-rho_sq *pow(x_scaled[i] - x_new_scaled[j],2));
  K2 = K' *Vm1;


  eps = mu0[1];
  sig = mu0[2];

  mu_M0 = mu_M + Kr2 * (y-mu_M);
  resid = y - mu_M0;

  neff = trace(Kr2);
  br = quad_form(Vm1,resid) / (N - 4 - neff);
  
  for (i in 1:N)
    y_pred_cont[i] = normal_rng(mu_M0[i],uy[i]);

  y_gp = K2 * (y-mu_M);
  y_m = phys_mod_vec(x_new,eps,sig,inadequacy); 
  y_conf = y_m + y_gp;

  for (i in 1:N2) {
    uyp = fmax(min(uy),y_conf[i] * mean(uy ./ y));
    y_pred[i] = normal_rng(y_conf[i],uyp);
  }
  

}


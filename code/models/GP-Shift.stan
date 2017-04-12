data {
  int N;
  vector[N] x;
  vector[N] y;
  vector[N] uy;
  int tag[N];
  int N2;
  vector[N2] x_new;
  int<lower=0, upper=1> inadequacy;
  vector[3] primu;
  cov_matrix[3] covprimu;
}
transformed data {
  int<lower=1> n_grp;  // nb groups
  real v_scale;
  vector[N] x_scaled;
  vector[N2] x_new_scaled;
  
  n_grp = max(tag);
  
  // Rescale x and variance for GP
  x_scaled = (x-min(x))/(max(x)-min(x));
  x_new_scaled = (x_new-min(x))/(max(x)-min(x));
  v_scale = max(uy)^2;
  
}
parameters {
  vector[3] mu0;
  real<lower= 0> eta_sq;
  real<lower= 0> inv_rho_sq;
  vector[n_grp] shift;
}
transformed parameters {
  real rho_sq;
  cov_matrix[N] V;
  vector[N] mu_M;
  vector[n_grp] shiftc;
  
  shiftc = shift - mean(shift);
  
  rho_sq = inv(inv_rho_sq);
  
  for (i in 1:N) 
    for (j in 1:N) 
      V[i,j] = eta_sq * v_scale * exp(-rho_sq * pow(x_scaled[i] - x_scaled[j],2)) 
             + (i==j ? uy[i]^2 : 0.0);
  
  for (i in 1:N) 
    mu_M[i] = phys_mod(x[i],mu0[1],mu0[2],inadequacy) 
            + shiftc[tag[i]];
  
}
model {
  
  mu0 ~ multi_normal(primu,covprimu);
  to_vector(shift) ~ normal(0,mu0[3]);
  
  eta_sq ~ normal(0,1);
  inv_rho_sq ~ normal(0,1);
  
  y ~ multi_normal(mu_M, V);

}
generated quantities{
  vector[N] resid;
  vector[N] mu_M0;
  vector[N] y_pred_cont;
  vector[N2] y_m;
  vector[N2] y_gp;
  vector[N2] y_conf;
  vector[N2] y_pred;
  real vy;
  real uyp;
  real eps;
  real sig;
  real tau;
  real neff;
  matrix[N,N]   Vm1;
  matrix[N,N]   Kr;
  matrix[N,N]   Kr2;
  matrix[N,N2]  K;
  matrix[N2,N]  K2;
  real br;

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
  tau = mu0[3];

  mu_M0 = mu_M + Kr2 * (y-mu_M);
  resid = y - mu_M0;

  neff = trace(Kr2);
  br = quad_form(Vm1,resid) / (N - 4 - neff - n_grp);
  
  for (i in 1:N)
    y_pred_cont[i] = normal_rng(mu_M0[i],uy[i]);

  y_gp = K2 * (y-mu_M);
  y_m = phys_mod_vec(x_new,eps,sig,inadequacy); 
  y_conf = y_m + y_gp;

  for (i in 1:N2) {
    uyp = fmax(min(uy),y_conf[i] * mean(uy ./ y));
    y_pred[i] = normal_rng(y_conf[i],(uyp^2+tau^2)^0.5);
  }

}


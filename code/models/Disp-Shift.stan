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
  real<lower= -4> log_sigma;
  real<lower= -4, upper=0> log_tau;
  vector[n_grp] shift;
}
transformed parameters {
  vector[N] mu_M;
  vector[N] V_T;
  vector[n_grp] shiftc;
  real sigma;
  real tau;
  
  shiftc = shift - mean(shift);
  
  sigma = 10^log_sigma;
  tau   = 10^log_tau;
  
  for (i in 1:N) { 
    mu_M[i] = phys_mod(x[i],eps,sig,inadequacy) 
            + shiftc[tag[i]];
    V_T[i] = uy[i]^2 + sigma^2;    
  } 
}
model {

  to_vector(shift) ~ normal(0,tau);
  
  for (i in 1:N) 
    y[i] ~ normal( mu_M[i],V_T[i]^0.5);
    
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
  
  br = quad_form(inverse(diag_matrix(V_T)),resid)  / (N-3-n_grp);
  
  for (i in 1:N)
    y_pred_cont[i] = normal_rng(mu_M[i],V_T[i]^0.5);

  for (i in 1:N2) {
    vy  = phys_mod(x_new[i],eps,sig,inadequacy);
    uyp = fmax(min(uy),vy * mean(uy ./ y));
    y_conf[i] = normal_rng(vy,sigma);
    y_pred[i] = normal_rng(y_conf[i], (uyp^2+tau^2)^0.5);
  }
}

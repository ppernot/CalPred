# Estimate prior for eps and sig from Disp results
load(paste0('simulation/',case,'/fitDisp.rda'))
eps = extract(fit,'eps')[[1]]
sig = extract(fit,'sig')[[1]]
primu = c(mean(eps),mean(sig))
cov_primu = cov(cbind(eps,sig))
rm(fit)

parOpt = c('eps', 'sig', 'eta_sq','inv_rho_sq')

# Provide reasonable starting point to shorten warmup
initf1 <- function() {
  list(mu0=primu, eta_sq = 0.5, inv_rho_sq = 0.5)
}

fit = stan(model_code = paste0(phys_model,stan_model), 
           model_name = model_tag,
           data = list(N =length(data$x), x=data$x, y=data$y, uy=data$uy,
                       N2=length(data$x_pred), x_new=data$x_pred,
                       inadequacy = inadequacy,
                       primu= primu, covprimu = cov_primu),
           pars = c(parOpt,'y_conf','y_pred','resid','br','y_gp','y_m',
                    'y_pred_cont','log_lik','neff'),
           init = initf1,
           iter = nb_iter, chains = nb_chains,
           # iter=2000, chains=1,
           warmup = nb_warmup, verbose=FALSE)

# library(loo)
# log_lik1 <- extract_log_lik(fit)
# loo1 <- loo(log_lik1)
# print(loo1,digits=3)
# print(mean(extract(fit,'neff')[[1]]))

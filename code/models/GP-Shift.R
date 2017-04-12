# Estimate prior for eps, sig and tau from Disp-Shift results
load(paste0('simulation/',case,'/fitDisp-Shift.rda'))
eps = extract(fit,'eps')[[1]]
sig = extract(fit,'sig')[[1]]
tau = extract(fit,'tau')[[1]]
primu = c(mean(eps),mean(sig),mean(tau))
cov_primu = cov(cbind(eps,sig,tau))
rm(fit)

parOpt = c('eps', 'sig', 'tau', 'eta_sq', 'inv_rho_sq')

# Provide reasonable starting point to shorten warmup
initf1 <- function() {
  list(mu0=primu, 
       eta_sq = 0.5, 
       inv_rho_sq = 0.5)
}

fit = stan(model_code = paste0(phys_model,stan_model), 
           model_name = model_tag,
           data = list(N =length(data$x), x=data$x, y=data$y, 
                       uy=data$uy, tag = data$tag,
                       N2=length(data$x_pred), x_new=data$x_pred,
                       inadequacy = inadequacy,
                       primu= primu, covprimu = cov_primu),
           pars = c(parOpt,'y_conf','y_pred','resid','br',
                    'y_gp','y_m','shiftc','y_pred_cont'),
           init = initf1,
           iter = nb_iter, chains = nb_chains,
           warmup = nb_warmup, verbose=FALSE)


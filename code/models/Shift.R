# Estimate prior for eps and sig from Disp results
load(paste0('simulation/',case,'/fitDisp.rda'))
eps = extract(fit,'eps')[[1]]
sig = extract(fit,'sig')[[1]]
tau = extract(fit,'sigma')[[1]]
rm(fit)

initf1 <- function() {
  list(eps   = mean(eps), 
       sig   = mean(sig), 
       log_tau  = log10(mean(tau)),
       shift = rep(0,max(data$tag))
       )
}
parOpt = c('eps', 'sig', 'tau')

fit = stan(model_code = paste0(phys_model,stan_model), 
           model_name = model_tag,
           data = list(N =length(data$x), x=data$x, y=data$y, uy=data$uy,
                       N2=length(data$x_pred), x_new=data$x_pred, tag = data$tag,
                       inadequacy = inadequacy),
           init = initf1,
           pars = c(parOpt,'y_conf','y_pred','resid','br','shiftc',
                    'y_pred_cont'),
           iter = nb_iter, chains = nb_chains, 
           warmup = nb_warmup, verbose=FALSE)



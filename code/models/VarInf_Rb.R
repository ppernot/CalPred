
# Estimate sigma from WLS results
load(paste0('simulation/',case,'/fitWLS.rda'))
eps = extract(fit,'eps')[[1]]
sig = extract(fit,'sig')[[1]]

lp  = extract(fit,'lp__')[[1]]
map = which.max(lp)
Rb  = extract(fit,'br')[[1]][map]
Ts  = (length(data$x)-2) / 2 * Rb # 2 parameters for WLS model
print(Ts)
rm(fit)

initf1 <- function() {
  list(eps   = mean(eps), 
       sig   = mean(sig)
  )
}
parOpt = c('eps', 'sig')

fit = stan(model_code = paste0(phys_model,stan_model), 
           model_name = model_tag,
           data = list(N =length(data$x), x=data$x, y=data$y, uy=data$uy,
                       N2=length(data$x_pred), x_new=data$x_pred,
                       inadequacy = inadequacy,
                       Ts = Ts),
           pars = c(parOpt,'y_conf','y_pred','resid','br',
                    'y_pred_cont'),
           init = initf1,
           iter = nb_iter, chains = nb_chains, 
           warmup = nb_warmup, verbose=FALSE)



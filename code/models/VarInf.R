
# Estimate sigma from Std results
load(paste0('simulation/',case,'/fitStd.rda'))
eps = extract(fit,'eps')[[1]]
sig = extract(fit,'sig')[[1]]
sigma0 = mean(extract(fit,'sigma')[[1]])

mse   = mean(apply(extract(fit,'resid')[[1]],2,mean)^2)
u_m2  = mean(apply(extract(fit,'y_conf')[[1]],2,sd)^2)
sigma = sigma0 * (mse / u_m2)^0.5

# lp = extract(fit,'lp__')[[1]]
# map = which.max(lp)
# Rb = extract(fit,'br')[[1]][map]
# T = (length(data$x)-3)/3*Rb # 3 parameters for Std model
# sigma1 = sigma0*T^0.5

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
                       sigma = sigma),
           pars = c(parOpt,'y_conf','y_pred','resid','br'),
           init = initf1,
           iter = nb_iter, chains = nb_chains, 
           warmup = nb_warmup, verbose=FALSE)



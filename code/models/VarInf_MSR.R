
# Estimate sigma from WLS results
load(paste0('simulation/',case,'/fitWLS.rda'))
eps = extract(fit,'eps')[[1]]
sig = extract(fit,'sig')[[1]]

lp  = extract(fit,'lp__')[[1]]
map = which.max(lp)

msr   = mean(extract(fit,'resid')[[1]][map,]^2)
u_m2  = mean(apply(extract(fit,'y_pred_cont')[[1]],2,sd)^2)
uy2   = mean(data$uy^2)
Ts    = (msr-uy2) / (u_m2-uy2)
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



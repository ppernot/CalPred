# Provide reasonable starting point to shorten warmup
load(paste0('simulation/',case,'/fitDisp.rda'))
eps = extract(fit,'eps')[[1]]
sig = extract(fit,'sig')[[1]]
rm(fit)

initf1 <- function() {
  list(eps = mean(eps), log_u_eps = 1+log10(sd(eps)), 
       sig = mean(sig), log_u_sig = 1+log10(sd(sig))
       )
}


parOpt = c('eps', 'u_eps', 'sig', 'u_sig','rho')

fit = stan(model_code = paste0(phys_model,stan_model), 
           model_name = model_tag,
           data = list(N =length(data$x), x=data$x, y=data$y, uy=data$uy,
                       N2=length(data$x_pred), x_new=data$x_pred,
                       inadequacy = inadequacy),
           pars = c(parOpt,'y_conf','y_pred','resid','br','mup',
                    'y_pred_cont'),
           control=list(adapt_delta=0.99, max_treedepth=10),
           init = initf1,
           iter = nb_iter, chains = nb_chains, 
           warmup = nb_warmup, verbose=FALSE)

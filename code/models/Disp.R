initf1 <- function() {
  list(eps   = 200, 
       sig   = 3.5, 
       log_sigma = -1
  )
}

parOpt = c('eps', 'sig', 'sigma')

fit = stan(model_code = paste0(phys_model,stan_model), 
           model_name = model_tag,
           data = list(N =length(data$x), x=data$x, y=data$y, uy=data$uy,
                       N2=length(data$x_pred), x_new=data$x_pred,
                       inadequacy = inadequacy, shift = shift),
           pars = c(parOpt,'y_conf','y_pred','resid','br','y_pred_cont'),
           init = initf1,
           iter = nb_iter, chains = nb_chains, 
           warmup = nb_warmup, verbose=FALSE)




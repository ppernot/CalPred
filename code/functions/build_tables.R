# Extract info from stanfit files to build tables

# Four types of samples with different treatments
# 1/ Get mean_over_sample and sd_over_sample
listInfo1=c('eps','sig','u_eps','u_sig','rho','sigma','tau')
# 2/ Get mean_over_control(mean_over_sample(.))
listInfo2=c('resid') # RMSD
# 3/ Get sqrt(mean_over_control(var_over_sample(.)))
listInfo3=c('y_pred_cont')
# 4/ Get value at MAP
listInfo4=c('br')
# 5/ Get mean_over_control(MAP)
listInfo5=c('resid') # MR

tab = data.frame()
for (model_tag in list_meth) {
  # Load data for current model
  fit_tag = paste0("fit",model_tag)
  load(paste0('simulation/',case,'/',fit_tag,'.rda'))
  
  # Process data lists
  sum=summary(fit)
  icol=0

  # 1/ Get mean_over_sample and sd_over_sample
  for (info in listInfo1) {
    icol=icol+1
    if(info %in% rownames(sum$summary))
      tab[model_tag,icol] = prettyUnc(sum$summary[info,'mean'],
                                     sum$summary[info,'sd'])
    else
      tab[model_tag,icol] =  '-'
  }
  
  # 2/ Get mean_over_control(mean_over_sample(.))
  for (info in listInfo2) {
    icol=icol+1
    X = try(extract(fit,info)[[1]],silent=TRUE)
    if(class(X) != 'try-error')
      tab[model_tag,icol] =  signif(mean(apply(X,2,mean)^2)^0.5,2)
    else
      tab[model_tag,icol] =  '-'
  }
  
  # 3/ Get sqrt(mean_over_control(var_over_sample(.)))
  for (info in listInfo3) {
    icol=icol+1
    X = try(extract(fit,info)[[1]],silent=TRUE)
    if(class(X) != 'try-error')
      tab[model_tag,icol] =  signif( mean( apply(X,2,var) )^0.5 ,2)
    else
      tab[model_tag,icol] =  '-'
  }
  
  # 4/ Get value at MAP
  for (info in listInfo4) {
    icol=icol+1
    X = try(extract(fit,info)[[1]],silent=TRUE)
    lp = extract(fit,'lp__')[[1]]
    map = which.max(lp)
    if(class(X) != 'try-error')
      tab[model_tag,icol] =  signif(X[map],2)
    else
      tab[model_tag,icol] =  '-'
  }
  
  # 5/ Get mean_over_control(MAP)
  for (info in listInfo5) {
    icol=icol+1
    X = try(extract(fit,info)[[1]],silent=TRUE)
    lp = extract(fit,'lp__')[[1]]
    map = which.max(lp)
    if(class(X) != 'try-error')
      tab[model_tag,icol] =  signif(mean(X[map,]),2)
    else
      tab[model_tag,icol] =  '-'
  }

  # # Mean time per efficient sample
  # icol = icol+1
  # tab[model_tag,icol] =  signif(sum(eltim[,2])/
  #                               nrow(X), 
  #                               5)  
  
}

# Tidy-up table, print and save
names(tab) = c(listInfo1,'RMSD','<u_e/D>','RB','MR')
orderCols = c(1:7,11,8,10,9)

print(knitr::kable(tab[,orderCols],'markdown',digits=2))

dir_out = paste0('simulation/',case,'/tables')
sink(file=paste0(dir_out,'/table_',case,'.tex'))
print(knitr::kable(tab[,orderCols],'latex',digits=2))
sink()

# Clean Environment
rm(X,tab)

print(fit, pars=c(parOpt,'br','lp__'))

sink(file=paste0('simulation/',case,'/diagnostics/summary_',model_tag,'.txt'))
print(fit, pars=c(parOpt,'br','lp__'))
cat('\n')
cat('*** Elapsed Time ***\n')
print(knitr::kable(eltim,format = 'markdown'))
sink()

tiff(file=paste0('simulation/',case,'/diagnostics/traces_',model_tag,'.tiff'),
    width=reso,height=reso)
print(traceplot(fit, pars=c(parOpt,'lp__'), 
                inc_warmup=TRUE, cex=2))
dev.off()

tiff(file=paste0('simulation/',case,'/diagnostics/pairs_',model_tag,'.tiff'),
    width=reso,height=reso)
pairs(fit, pars=c(parOpt,'lp__'), gap=0, cex=2)
dev.off()
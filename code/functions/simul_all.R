# Create arborescence (if absent) ####
dir.create(path='simulation', 
           showWarnings = FALSE)
dir.create(path=paste0('simulation/',case), 
           showWarnings = FALSE)
dir.create(path=paste0('simulation/',case,'/tables'), 
           showWarnings = FALSE)
dir.create(path=paste0('simulation/',case,'/figures'), 
           showWarnings = FALSE)
dir.create(path=paste0('simulation/',case,'/diagnostics'), 
           showWarnings = FALSE)

# Generate/load data ####
set.seed(repeat_seed)
if(substr(case,1,2)=='SD') {
  data = gen_synth_data(shift=shift)
} else {
  data = get_ref_data()
}

# Save a copy in latex format
sink(paste0('simulation/',case,'/data.tex'))
print(knitr::kable(cbind(data$tag,data$x,data$y,signif(data$uy,2)), 
                   format='latex'))
sink()

# Run simulations ####
for( model_tag in list_meth ) {
  simul_file = paste0('simulation/',case,'/fit',model_tag,'.rda')
  if(!file.exists(simul_file)) { 
    # Do not run if simul_file exists already
    stan_model = 
      paste0(
        readLines(con=paste0('code/models/',model_tag,'.stan')),
        collapse='\n'
      )
    source(paste0('code/models/',model_tag,'.R'))
    eltim = get_elapsed_time(fit)
    
    # Save fit object and context
    save(fit, data, parOpt, eltim, file=simul_file)

    # Generate diagnostic outputs
    source('code/functions/diagnostics.R')
    
    # Clean environment of big object
    rm(fit)
  }
}

# Generate figures ####
source('code/functions/build_figures.R')

# Generate table(s) ####
source('code/functions/build_tables.R')

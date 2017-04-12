# Extract info from stanfit files to build figures

for (model_tag in list_meth) {
  # Load data for current model
  fit_tag = paste0("fit",model_tag)
  load(paste0('simulation/',case,'/',fit_tag,'.rda'))
  
  # Plot
  plot_all(legend=list_legends[model_tag])
  
}
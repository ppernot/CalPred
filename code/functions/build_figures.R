# Extract info from stanfit files to build figures

for (model_tag in list_meth) {
  # Load data for current model
  fit_tag = paste0("fit",model_tag)
  load(paste0('simulation/',case,'/',fit_tag,'.rda'))
  
  # Plot
  plot_all(legend=list_legends[model_tag])
  
}

## Plot shifts for relevant methods

list_meth_shift = list_meth[grep('Shift',list_meth)]

if(length(list_meth_shift) > 0)
  plot_shifts(model_tags=list_meth_shift, 
              true_shifts=data$true_shifts, 
              cexPlot = cexPlot)
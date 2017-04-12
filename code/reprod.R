# Clean environment
rm(list = ls())

# Set misc. options and load libraries and functions ####
options(width=120)

libs = c('parallel', 'rstan', 'knitr', 'aplpack', 'RColorBrewer')
void = sapply(libs,
              function(x)
                if(!require(x,
                            character.only = T,
                            warn.conflicts = F,
                            quietly = T)
                ) 
                  install.packages(x)
)
rm(void,libs)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Physical model functions for stan 
phys_model = 
  paste0(
    readLines(con='code/models/phys_model.stan'),
    collapse='\n'
  )

# Misc. R functions
source('code/models/phys_model.R')
source('code/functions/gen_synth_data.R')
source('code/functions/plot_funcs.R')

# Sampling options
nb_warmup = 1000
nb_iter   = 5000
nb_chains = 4

#################################################################

# Case: SD-1 ####################################################
case= 'SD-1'
inadequacy = TRUE  # Model inadequacy
shift      = FALSE # Data inconsistency

## Define methods
list_meth=c('Std','Disp','GP','WLS','VarInf_Rb',
            'VarInf_MSR','Margin','ABC','HierC','HierC-loc')[1:10]
list_legends=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)
names(list_legends)= list_meth

## Define plot limits
eps_lim = c(180,300)     # Limits for epsilon plots
sig_lim = c(3.45, 3.65)  # Limits for sigma plots
res_lim = c(-1,1) * 0.6  # Limits for residuals plots

## Run simulations 
source('code/functions/simul_all.R')

# Case: SD-2 ####################################################
case= 'SD-2'
inadequacy = FALSE
shift      = TRUE

## Define methods
list_meth=c('Disp','Shift','Cov','Hier','Hier-Shift','Hier-Cov')
list_legends=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE)
names(list_legends)= list_meth

## Define plot limits 
eps_lim = c(170,230)     # Limits for epsilon plots
sig_lim = c(3.58,3.62)   # Limits for sigma plots
res_lim = c(-1,1) * 1.5  # Limits for residuals plots

## Run simulations for the current case
source('code/functions/simul_all.R')

## Plot shifts for relevant methods
plot_shifts(model_tags=list_meth[c(2,5)], 
            true_shifts=data$true_shifts, 
            cexPlot = cexPlot)

# Case: SD-3 ####################################################
case= 'SD-3'
inadequacy = TRUE
shift      = TRUE

## Define methods
list_meth=c('Disp','Disp-Shift','GP-Shift','Margin-Shift',
            'ABC-Shift','Hier-Shift','Hier-Cov')
list_legends=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
names(list_legends)= list_meth

## Define plot limits 
eps_lim = c(180,280)     # Limits for epsilon plots
sig_lim = c(3.50,3.62)   # Limits for sigma plots
res_lim = c(-1,1) * 1.0  # Limits for residuals plots

## Run simulations for the current case
source('code/functions/simul_all.R')

## Plot shifts for relevant methods
plot_shifts(model_tags=list_meth[grep('Shift',list_meth)], 
            true_shifts=data$true_shifts, 
            cexPlot = cexPlot)

# Case: Kr ######################################################
case= 'Kr'
inadequacy = FALSE
shift      = TRUE

## Define methods
list_meth=c('Disp','Disp-Shift','GP-Shift','Margin-Shift',
            'ABC-Shift')
list_legends=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
names(list_legends)= list_meth

## Define plot limits
eps_lim = c(180,210)     # Limits for epsilon plots
sig_lim = c(3.52,3.58)   # Limits for sigma plots
res_lim = c(-1,1) * 1.5  # Limits for residuals plots

## Run simulations for the current case
source('code/functions/simul_all.R')

## Plot shifts for relevant methods
plot_shifts(model_tags=list_meth[grep('Shift',list_meth)], 
            cexPlot = cexPlot)

# Session Info ####
sink(file='session_info.txt')
writeLines(readLines(file.path(Sys.getenv("HOME"), ".R/Makevars")))
devtools::session_info()
sink()

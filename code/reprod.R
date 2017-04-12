# Clean environment
rm(list = ls())

# Set misc. options and load libraries and functions ####
options(width=120)

if(!require(knitr))
  install.packages("knitr",dependencies=TRUE)
library(knitr) # Provides kable

if(!require(parallel))
  install.packages("parallel",dependencies=TRUE)
library(parallel)
options(mc.cores = parallel::detectCores())

if(!require(rstan))
  install.packages("rstan",dependencies=TRUE)
library(rstan)
rstan_options(auto_write = TRUE)

# Physical model functions for stan and R 
phys_model = 
  paste0(
    readLines(con='code/models/phys_model.stan'),
    collapse='\n'
  )
source('code/models/phys_model.R')

# Misc. R functions
source('code/functions/gen_synth_data.R')
source('code/functions/plot_funcs.R')

# Sampling parameters
nb_warmup   = 1000
nb_iter     = 5000
nb_chains   = 4
repeat_seed = 127   # Random seed used for repeatability

#################################################################

# Case: SD-1 ####################################################
case= 'SD-1'
inadequacy = TRUE  # Model inadequacy
shift      = FALSE # Data inconsistency

## Define methods
list_meth=c('Std','Disp','GP','WLS','VarInf_Rb','VarInf_MSR','Margin','ABC','HierC')
list_legends=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)
names(list_legends)= list_meth

## Define plot limits
eps_lim = c(180,300)     # Limits for epsilon plots
sig_lim = c(3.45, 3.65)   # Limits for sigma plots
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

# Session Info ####
sink(file='session_info.txt')
print(sessionInfo())
sink()

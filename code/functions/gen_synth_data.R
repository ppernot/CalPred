Ferror = function(x,F0,g) {
  # Model of heteroscedastic measurement uncertainty
  return( F0 * exp( g * abs(1/x-1/300) ) )
}

gen_synth_data = function (shift = FALSE) {
  # Reference values of parameters
  eps0 = 195
  sig0 = 3.6
  true_point = c(eps0,sig0)
  
  # Tempetature range
  xMin = 120
  xMax = 2000
  
  # Uncertainty range (rel.)
  F0 = 0.001
  g  = 50

  # Syst. error sd
  if (shift)
    uShift = 0.5
  else
    uShift = 0
  
  # Number of series
  nSeries = 10
  
  # Series sizes range
  nMin = 5 ; nMax= 15
  
  ipt = 0
  x = y = uy = err = tag = true_shifts = c()
  for (iSeries in 1:nSeries) {

    # Size of series
    n = sample(nMin:nMax,1)
    
    # Select temp range for series
    x0 = runif(1, min=xMin, max=xMax/2)
    x1 = runif(1, min=x0, max=xMax)
    
    # Systematic error
    shift = uShift * rnorm(1, mean=0, sd=1)
    true_shifts[iSeries] = shift
    
    # Step size in series
    dx = (x1-x0)/(n-1)
    for (j in 1:n) {
      ipt = ipt +1
      x[ipt]  = x0 + (j-1)*dx

      y0 = phys_mod(x[ipt],eps0,sig0, inadequacy=FALSE)
      
      uy[ipt] = y0 * Ferror(x[ipt],F0,g)
      y[ipt]  = y0 + shift +  rnorm(1, mean=0, sd=uy[ipt])
      err[ipt]= y[ipt] - y0
      tag[ipt] = ifelse(shift,iSeries,1)  
    }
  }

  # Define prediction grid
  x_pred=seq(min(x)*0.5,max(x)*1.1,len=100)   
  
  return(
    list(
      x=x, y=y, uy=uy,
      tag=tag, err=err,
      true_point=true_point,
      true_shifts=true_shifts,
      x_pred = x_pred
    )
  )
}

get_ref_data = function() {
  switch(case,
         Kr = { dc  = read.csv(file='data/data_Kr_visc_Bich1990.csv',
                               comment.char = "#")
                x   = dc[,1]
                y   = dc[,2]
                uy  = dc[,3] / 100 * y # Convert to absolute uncert
                tag = dc[,4]
                cel = dc[,5]
                x = x + cel*273.15 # Convert to Kelvin if necessary
              },
         Ar = { dc = read.table(file='data/data_Ar_visc_pri.csv', 
                                stringsAsFactors = FALSE)
                x   = dc[,1]
                y   = dc[,2]/10
                uy  = dc[,3]/10
                tag = dc[,4]
                tag=tag-min(tag)+1
         }
  )
  
  # Define prediction grid
  x_pred=seq(min(x)*0.5,max(x)*1.1,len=100)   
  
  return(
    list(
      x=x, y=y, uy=uy,
      tag=tag, err=NULL,
      true_point=NULL,
      true_shifts=NULL,
      x_pred = x_pred
    )
  )  
}
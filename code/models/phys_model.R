phys_mod = function(Temp,eps,sig,inadequacy=FALSE) {
  # Chapman-Enskog viscosity model
  A = 1.16145; B =-0.14874; C = 0.52487
  D =-0.77320; E = 2.16178; F =-2.43787
  
  # Introduce inadequacy
  if(inadequacy) C = C*0.5
  
  # M = 39.948 # Ar (g/mol)
  M = 83.798 # Kr (g/mol)
  Tstar = Temp/eps
  omega = A*Tstar^B + C*exp(D*Tstar) + E*exp(F*Tstar)
  # Output in microPa.s
  return( 2.6693 * (M*Temp)^0.5 / (sig^2 * omega) )
}



functions {
  real phys_mod(real x, real eps, real sig, int inadequacy) {
    real Tstar;
    real omega;
    real A; 
    real B; 
    real C; 
    real D;
    real E;
    real F;
    real M;
    real eps0;
    real sig0;
    
    A =  1.16145;
    B = -0.14874;
    C =  0.52487;
    D = -0.77320;
    E =  2.16178;
    F = -2.43787;
    M =  83.798;
    
    if (inadequacy == 1)
      C = C*0.5;
      
    eps0 = fmax(eps,0.1);
    sig0 = fmax(sig,0.01);
    
    Tstar = x / eps0;
    omega = A*exp(B*log(Tstar)) + C*exp(D*Tstar) + E*exp(F*Tstar);
    return 2.6693 * exp(0.5*log(M*x)) / (sig0^2 * omega);
  }
  vector phys_mod_vec(vector x, real eps, real sig, int inadequacy) {
    vector[rows(x)] tmp;
    for (i in 1:rows(x))
      tmp[i] = phys_mod(x[i],eps,sig,inadequacy);
    return tmp;
  }
  vector grad_phys_mod(real x, real eps, real sig, int inadequacy) {
    vector[2] grad;
    real A; 
    real B; 
    real C; 
    real D;
    real E;
    real F;
    real M;
    real eps0;
    real sig0;

    real expr3;
    real expr4;
    real expr5;
    real expr9;
    real expr13;
    real expr15;
    real expr16;
    real expr19;
    real expr35;

    A =  1.16145;
    B = -0.14874;
    C =  0.52487;
    D = -0.77320;
    E =  2.16178;
    F = -2.43787;
    M =  83.798;

    if (inadequacy == 1)
      C = C*0.5;

    eps0 = fmax(eps,0.1);
    sig0 = fmax(sig,0.01);
    
    expr3 = 2.6693 * (M * x)^0.5;
    expr4 = sig0^2;
    expr5 = x/eps0;
    expr9 = exp(D * expr5);
    expr13 = exp(F * expr5);
    expr15 = A * expr5^B + C * expr9 + E * expr13;
    expr16 = expr4 * expr15;
    expr19 = x/eps0^2;
    expr35 = expr16^2;
  
    grad[1] = expr3 * 
                (expr4 * 
                 (E * (expr13 * (F * expr19)) + 
                       (C * (expr9 * (D * expr19)) + 
                        A * (expr5^(B - 1) * (B * expr19)))))/expr35;
    grad[2] = -(expr3 * (2 * sig0 * expr15)/expr35);
    return grad;
  }
  matrix grad_phys_mod_vec(vector x, real eps, real sig, int inadequacy) {
    matrix[rows(x),2] tmp;
    for (i in 1:rows(x))
      tmp[i,1:2] = grad_phys_mod(x[i],eps,sig,inadequacy)';
    return tmp;
  }
}

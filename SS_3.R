# Calculate steady state for model 3 - see appendix to paper for derivation of this

# Constants used in the solution
Y <- (tau+u+m)

if (q>0){ # With reinfection
  # Parts of quadratic in I
  cq2 <- u*u*(Y-a*beta) + c*u*(Y-beta-tau)
  bq2 <- Y*beta*(q*a*u+c+u) - a*beta*u*(m+q*(beta + tau)) - c*beta*(m+tau)
  aq2 <- a*q*beta*beta*(Y-m-tau)
  # Solution to quadratic for I
  I <- max(c(0,(-bq2+sqrt(bq2*bq2 - 4*aq2*cq2))/(2*aq2),(-bq2-sqrt(bq2*bq2 - 4*aq2*cq2))/(2*aq2)))
  # Steady state solution for S and L
  S <- (u+m*I)/(u+beta*I)
  L <- ((1-a)*beta*I*S+tau*I)/(a*q*beta*I+c+u)
  
  if (a==0){ # if a=0 (as in ragonnet parameterisation) then reinfection has no effect, the steady state calculation also needs redoing as aq2 = 0 when a =0
    I <- max(0,(c*beta*u - u*(Y*(c+u)-c*tau))/((Y*(c+u)-c*tau)*beta - c*beta*m))
    S <- (u+m*I)/(u+beta*I)
    L <- (beta*I*S+tau*I)/(c+u)
  }
  
}
if (q==0){ # No reinfection
  I <- max(0,((c+u)*u*Y - (c+u)*a*beta*u - c*(1-a)*beta*u - c*tau*u)/(c*(1-a)*beta*m + c*tau*beta - (c+u)*beta*Y + (c+u)*a*beta*m))
  S <- (u+m*I)/(u+beta*I)
  L <- ((1-a)*(beta*S*I)+tau*I)/(c+u) 
}
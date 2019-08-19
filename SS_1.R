# Calculate the steady state for model 1 - see appendix to paper for derivation of this

# Constants used in the solution
Z <- k+e+u
Y <- tau+u+m

if (q>0){ # with reinfection
  # Parts of quadratic in I
  cq2 <- u*Z*((Z*Y-beta*k)*(c+u) - c*(e*beta+Z*tau))
  bq2 <- (beta*Z*((Z*Y-m*k)*(c+u)-c*(e*m+tau*Z)) + u*beta*q*((Z*Y-beta*k)*(Z-e)-k*(e*beta+Z*tau)))
  aq2 <- q*beta*beta*((Z*Y-m*k)*(Z-e) - k*(e*m+tau*Z))
  # Solution to quadratic for I
  I <- max(c(0,(-bq2+sqrt(bq2*bq2 - 4*aq2*cq2))/(2*aq2),(-bq2-sqrt(bq2*bq2 - 4*aq2*cq2))/(2*aq2)))
  S <- (u+m*I)/(u+beta*I)
  LB <- (e*beta*I*S + Z*tau*I)/(Z*(q*beta*I+c+u) - e*q*beta*I)
  LA <- beta*I*(S+q*LB)/Z
}

if(q==0){ # no reinfection
  I <- max(0,(beta*k*u*(c+u) + c*beta*e*u + c*u*Z*tau - u*Z*Y*(c+u))/(beta*Z*Y*(c+u) - beta*k*m*(c+u) - c*beta*e*m - c*beta*Z*tau))
  S <- (u+m*I)/(u+beta*I)
  LA <- (beta*S*I)/Z 
  LB <- (e*LA + tau*I)/(c+u)
} 


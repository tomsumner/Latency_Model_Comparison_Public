# Calculate the steady state for model 2 - see appendix to paper for derivation of this

# Constants used in the solution
Y <- (tau+u+m)

if (q>0){ # with reinfection
  # Parts of quadratic in I
  cq2 <- beta*u*(k*(b*u+c)+c*u*(1-b))+(k+u)*u*(c*tau - Y*(c+u))
  bq2 <- k*beta*b*(m*(c + u) + q*u*(tau + beta)) + beta*(k+u)*(c*(m*(1-b)+tau) - Y*(q*b*u+c+u))
  aq2 <- -1*q*beta*beta*b*u*(Y+k)
  # Solution to quadratic for I
  I <- max(c(0,(-bq2+sqrt(bq2*bq2 - 4*aq2*cq2))/(2*aq2),(-bq2-sqrt(bq2*bq2 - 4*aq2*cq2))/(2*aq2)))
  S <- (u+m*I)/(u+beta*I)
  LB <- ((1-b)*beta*I*S+tau*I)/(q*beta*b*I+c+u)
  LA <- beta*I*b*(S+q*LB)/(k+u)
}
if (q==0){ # no reinfection
  I <- max(0,u*(Y*(k+u)*(c+u)-(k+u)*c*tau-(c+u)*k*b*beta-(k+u)*c*(1-b)*beta)/
             (m*((c+u)*k*b*beta+(k+u)*c*(1-b)*beta)-beta*(Y*(k+u)*(c+u)-(k+u)*c*tau)))
  S <- (u+m*I)/(u+beta*I)
  LA <- b*beta*S*I/(k+u) 
  LB <- ((1-b)*beta*I*S + tau*I)/(c+u)
}


# Model 2 

Cum_2 <- function (pars,state,times=seq(0,t_max)) {
  derivs <- function(t,state,pars){
    with(as.list(c(state,pars)),{
      
      # Latent early
      dLA <- -k*LA
      
      # Latent late
      dLB <- -c*LB
      
      # Infectious
      dI <- k*LA + c*LB                 
      
      # Return these things
      return(list(c(dLA, dLB, dI)))
    })
  }
  
  ## ode solves the model by integration...
  return(ode(y = state, times = times, func = derivs, parms = pars))
}
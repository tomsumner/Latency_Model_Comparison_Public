# Model 1 

Cum_1 <- function (pars,state,times=seq(0,t_max)) {
  derivs <- function(t,state,pars){
    with(as.list(c(state,pars)),{
    
      # Latent early
      dLA <- - k*LA - e*LA

      # Latent late
      dLB <- e*LA - c*LB
    
      # Infectious
      dI <- k*LA + c*LB                

      # Return these things
      return(list(c(dLA, dLB, dI)))
    })
  }
  
  ## ode solves the model by integration...
  return(ode(y = state, times = times, func = derivs, parms = pars))
}
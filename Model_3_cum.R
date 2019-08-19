# Model 3 

Cum_3 <- function (pars,state,times=seq(0,t_max)) {
  derivs <- function(t,state,pars){
    with(as.list(c(state,pars)),{
      
      # Latent 
      dL <- -c*L 
      # Infectious
      dI <- c*L 
      
      # Return these things
      return(list(c(dL, dI)))
    })
  }
  
  ## ode solves the model by integration...
  return(ode(y = state, times = times, func = derivs, parms = pars))
}
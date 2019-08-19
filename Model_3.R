# Model 3 with intervention. 

TB_model_3 <- function (pars,state,times=seq(0,tend)) {
  derivs <- function(t,state,pars){
    with(as.list(c(state,pars)),{
      
      # detection rate - allows it to be changed to alter shape of incidence curve
      cov <- approxfun(x=c(1,2),y=c(0,cov_PT),method="linear",rule=2)  
  
      # Check total population stays at 1
      Total = S+L+I+P
      
      # Births set to total deaths to keep total population constant 
      births <- u*Total + m*I 
      
      # force of infection
      FOI <- beta*I
    
      # Total number given PT in this year
      Number_PT <- cov(t)*L
      
      # Susceptible
      dS <- births - u*S -                   # plus births minus deaths 
            FOI*S                            # minus infection
      
      # Latent 
      dL <- (1-a)*S*FOI -                    # plus infection S
            c*L -                            # minus progression from L to disease
            u*L -                            # minus death
            q*FOI*a*L +                      # minus reinfection in L
            (1-a)*q*FOI*P +                  # plus reinfection from P
            tau*I -                          # plus recovery from I
            cov(t)*L                         # minus those given effective PT
      
      # Infectious
      dI <- FOI*S*a +                        # primary disease
            c*L +                            # reactivation from L
            w*c*P +                          # progression from P
            FOI*q*a*(L+P) -                  # exogenous reinfection
            tau*I - u*I - m*I                # minus recovery and death              
      
      # Post PT
      dP <- cov(t)*L -                       # plus got PT
            w*c*P -                          # minus progression to disease - if W=0, then PT is completley protective against reactivation
            q*P*FOI -                        # minus reinfection
            u*P                              # minus death
        
      # keep track of cumulative number given PT
      dN_PT <- Number_PT
      
      # Derived outputs
      Inc <- c*L + w*c*P + FOI*S*a + FOI*q*a*(L+P)   # Total Incidence
      
      # Return these things
      return(list(c(dS, dL, dI, dP, dN_PT), 
                  Total = Total,
                  Inc = Inc,
                  Number_PT = Number_PT))
    })
  }
  
  ## ode solves the model by integration...
  return(ode(y = state, times = times, func = derivs, parms = pars))
}
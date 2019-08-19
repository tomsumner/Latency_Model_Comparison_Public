# Model 2 with intervention. Includes a single post PT state

TB_model_2a <- function (pars,state,times=seq(0,tend)) {
  derivs <- function(t,state,pars){
    with(as.list(c(state,pars)),{
      
      # detection rate - allows it to be changed to alter shape of incidence curve
      cov <- approxfun(x=c(1,2),y=c(0,cov_PT),method="linear",rule=2)  
  
      # Check total population stays at 1
      Total = S+LA+LB+I+P
      
      # Births set to total deaths to keep total population constant 
      births <- u*Total + m*I 
      
      # force of infection
      FOI <- beta*I
    
      # Total number given PT in this year
      Number_PT <- cov(t)*(LA+LB)
      
      # Susceptible
      dS <- births - u*S -                    # plus births minus deaths 
            FOI*S                             # minus infection
      
      # Latent early
      dLA <- b*S*FOI +                        # plus infection S
             q*b*(LB+P)*FOI -                 # plus reinfection
             k*LA -                           # minus progression from LA to disease
             u*LA -                           # minus death
             cov(t)*LA                        # minus those given effective PT
      
      # Latent late
      dLB <- (1-b)*S*FOI +                    # plus infection from S
             q*(1-b)*FOI*P -                  # plus reinfection
             c*LB -                           # minus progression from LB to disease    
             q*b*LB*FOI -                     # minus reinfection 
             u*LB -                           # minus death  
             cov(t)*LB +                      # minus those given effective PT
             tau*I                            # plus recovery from I                    
      
      # Infectious
      dI <- k*LA + c*LB +                     # plus progression from LA and LB
            w*c*P -                           # plus progression from PA and PB
            tau*I - u*I - m*I                 # minus recovery and death              
      
      # Post PT (early) 
      dP <- cov(t)*(LA+LB) -                  # plus got PT
            w*c*P -                           # minus progression to disease
            q*P*FOI -                         # minus reinfection
            u*P                               # minus death
      
      # keep track of cumulative number given PT
      dN_PT <- Number_PT
      
      # Derived outputs
      Inc <- k*LA + c*LB + w*c*P              # Total Incidence
      
      # Return these things
      return(list(c(dS, dLA, dLB, dI, dP, dN_PT), 
                  Total = Total,
                  Inc = Inc,
                  Number_PT = Number_PT))
    })
  }
  
  ## ode solves the model by integration...
  return(ode(y = state, times = times, func = derivs, parms = pars))
}
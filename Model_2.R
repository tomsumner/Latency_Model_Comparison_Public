# Model 2 with intervention. Includes 2 post PT states. 

TB_model_2 <- function (pars,state,times=seq(0,tend)) {
  derivs <- function(t,state,pars){
    with(as.list(c(state,pars)),{
      
      # coverage of prevetnive therpay - implemented as an approximating function
      cov <- approxfun(x=c(1,2),y=c(0,cov_PT),method="linear",rule=2)  
  
      # Check total population stays at 1
      Total = S+LA+LB+I+PA+PB
      
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
             q*b*(LB+PB)*FOI +                # plus reinfection
             q*PA*FOI -
             k*LA -                           # minus progression from LA to disease
             u*LA -                           # minus death
             cov(t)*LA                        # minus those given effective PT
      
      # Latent late
      dLB <- (1-b)*S*FOI +                    # plus infection from S
             q*(1-b)*FOI*PB -                 # plus reinfection
             c*LB -                           # minus progression from LB to disease    
             q*b*LB*FOI -                     # minus reinfection 
             u*LB -                           # minus death  
             cov(t)*LB +                      # minus those given effective PT
             tau*I                            # plus recovery from I                    
      
      # Infectious
      dI <- k*LA + c*LB +                     # plus progression from LA and LB
            w*(k*PA + c*PB) -                 # plus progression from PA and PB
            tau*I - u*I - m*I                 # minus recovery and death              
      
      # Post PT (early) 
      dPA <- cov(t)*LA -                      # plus got PT
             w*k*PA -                         # minus progression to disease
             q*PA*FOI -                       # minus reinfection
             u*PA                             # minus death
      
      # Post PT (early) 
      dPB<- cov(t)*LB -                       # plus got PT
            w*c*PB -                          # minus progression to disease
            q*PB*FOI -                        # minus reinfection
            u*PB                              # minus death
      
      # keep track of cumulative number given PT
      dN_PT <- Number_PT
      
      # Derived outputs
      Inc <- k*LA + c*LB + w*(k*PA + c*PB)    # Total Incidence
      
      # Return these things
      return(list(c(dS, dLA, dLB, dI, dPA, dPB, dN_PT), 
                  Total = Total,
                  Inc = Inc,
                  Number_PT = Number_PT))
    })
  }
  
  ## ode solves the model by integration...
  return(ode(y = state, times = times, func = derivs, parms = pars))
}
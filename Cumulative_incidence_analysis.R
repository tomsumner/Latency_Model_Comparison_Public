# CUMULATIVE_INCIDENCE_ANALYSIS.R

# This script compares the outputs from 3 different model structures using the parameters reported in Menzies and Ragonnet
# It looks at the cumulative TB incidence over time following infection

#######################################################################################################################
# DEFINE PARAMETERS

# We consider 2 parameter sets, those estimated by Menzies and those estimated by Ragonnet
# Ragonnet parameters are reported in daily units 

# These parameters are the same for all models
tau <- 1         # rate of recovery
u <- 1/50        # natural mortality
m <- 0.03        # TB mortality
qr <- c(0.5,0)   # rel sus after prior infection - consider 2 situations: 0=no reinfection; 0.5 50% protection against disease following reinfection

# These parameters vary by model and source (Menzies or Ragonnet)
# proportion direct TB - NOTE THE VALUE FOR MODEL 3 IN RAGONNET IS NOT REPORTED IN THE PAPER SO WAS OBTAINED BY FITTING
ar <- cbind(c(0,0,0.0665),c(0,0,0.085))  
# proportion to early latent state 
br <- cbind(c(0,0.0860,0),c(0,1-0.91,0))  
# progression from late latent state - NOTE THE VALUE FOR MODEL IN RAGONNET IS NOT REPORTED IN THE PAPER SO WAS OBTAINED BY FITTING
cr <- cbind(c(0.000594,0.000594,0.00337),c(5.5e-6,5.5e-6,8.52e-6))       
# progression from early latent state 
kr <- cbind(c(0.0826,0.955,0),c(1.1e-3,1.1e-2,0))   
# movement to late latent state 
er <- cbind(c(0.872,0,0),c(1.1e-2,0,0))                        

#######################################################################################################################
# CALCUATE AND PLOT CUMULATIVE INCIDENCE

pars <- c()   # no parameter to pass
t_max <- 7300 # this is 20 years in days

# Array to hold outputs
Cum_out_S <- array(NA,dim=c(t_max+1,dim(er)[1],dim(er)[2]))

for (ii in 1:dim(er)[1]){         # For each model 
  for (zz in 1:dim(er)[2]){       # For each parameterisation
    
    # Pick the right parameters
    a <- ar[ii,zz]
    b <- br[ii,zz]
    c <- cr[ii,zz]
    k <- kr[ii,zz]
    e <- er[ii,zz]
    
    # define the initial conditions
    state_cum_S <- list(c(LA=1,LB=0,I=0),                      # Model 1 - all in early latent state
                        c(LA=br[2,zz],LB=1-br[2,zz],I=0),  # Model 2 - split between early and late latent states
                        c(L=1-ar[3,zz],I=ar[3,zz]))        # Model 3 - split between latent and disease
    state <- unlist(state_cum_S[ii])
    
    # Run the model
    if (ii==1) temp <- as.data.frame(Cum_1(pars,state))
    if (ii==2) temp <- as.data.frame(Cum_2(pars,state))
    if (ii==3) temp <- as.data.frame(Cum_3(pars,state))
    
    # Store the outputs
    Cum_out_S[,ii,zz] <- temp$I
    
  }
  
}

# Plot results for both parameterisations
Cum_inc <- as.data.frame(rbind(cbind(Cum_out_S[,,1],seq(0,t_max),"A"),
                                   cbind(Cum_out_S[,,2],seq(0,t_max/365,by=1/365),"B"))) # This converts daily to annual
colnames(Cum_inc) <- c("1","2","3","Year","Source")
c_m <- melt(Cum_inc,id.vars=c("Year","Source"))

# Figure 1 in appendix 
Figure_A1 <- ggplot(data=c_m,aes(as.numeric(as.character(Year)),as.numeric(as.character(value))))+
            geom_line(data=c_m[c_m$Source%in%c("A","B"),],aes(color=variable,linetype=as.factor(Source)),size=1)+
  theme_bw()+ theme(legend.position = "bottom")+
  ylab(c("Cumulative incidence of TB"))+xlab(c("Years since infection"))+
  scale_colour_manual(name="",values=cbPalette)+
  scale_linetype_manual(values = c("dashed","solid"),name="")+
  scale_y_continuous(limits=c(0,0.15),expand=c(0,0))+
  scale_x_continuous(limits=c(0,20),expand=c(0,0))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  guides(col = guide_legend(nrow = 1))+
  guides(linetype = guide_legend(nrow = 1))

#######################################################################################################################
# CALCULATE TB RISKS FOR ALL PARAMETER OPTIONS

# Ragonnet parameters now converted to annual units       

# proportion direct TB
ar_S <- cbind(c(0,0,0.0665),c(0,0,0.085))    
# proportion to early latent state 
br_S <- cbind(c(0,0.0860,0),c(0,1-0.91,0))   
# progression from late latent state
cr_S <- cbind(c(0.000594,0.000594,0.00337),365*c(5.5e-6,5.5e-6,8.52e-6))
# progression from early latent state 
kr_S <- cbind(c(0.0826,0.955,0),365*c(1.1e-3,1.1e-2,0))
# movement to late latent state 
er_S <- cbind(c(0.872,0,0),365*c(1.1e-2,0,0)) 
                                                                  
# Life time risk
L_S <- rbind(
  # Model 1 
  (kr_S[1,]*(cr_S[1,]+u) + (er_S[1,]*cr_S[1,]))/((kr_S[1,]+er_S[1,]+u)*(cr_S[1,]+u)),
  # Model 2 
  (br_S[2,]*kr_S[2,]/(kr_S[2,]+u)) + (1-br_S[2,])*cr_S[2,]/(cr_S[2,]+u),
  # Model 3
  ar_S[3,] + (1-ar_S[3,])*(cr_S[3,]/(cr_S[3,]+u)))

# Risk in "early" latency
LE_S <- rbind(
  # Model 1
  kr_S[1,]/(kr_S[1,]+er_S[1,]+u),
  # Model 2
  (br_S[2,]*kr_S[2,]/(kr_S[2,]+u)),
  # Model 3
  ar_S[3,])

# Risk in "late" latency
LL_S <- rbind(
  # Model 1
  (er_S[1,]/(er_S[1,]+kr_S[1,]+u))*(cr_S[1,]/(cr_S[1,]+u)),
  # Model 2
  (1-br_S[2,])*cr_S[2,]/(cr_S[2,]+u),
  # Model 3
  (1-ar_S[3,])*cr_S[3,]/(cr_S[3,]+u))

# Proportion of risk in early latency 
per_early_S <- LE_S/L_S





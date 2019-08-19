# Analysis of model 1 with different duration of early latency

duration_r <- seq(1,5,by=1)      # duration of early latency 
er_duration <- 1/duration_r      # rate of move from early to late latency - 1/duration
qr_duration <- c(0.5,0)          # rel sus after prior infection - consider 2 situations: 0=no reinfection; 0.5 50% protection
c <- 0.000594                    # rate of progression from late latency - fixed to value from Menzies
tau <- 1                         # rate of recovery
u <- 1/50                        # natural mortality
m <- 0.03                        # TB mortality

# Set the lifetime risk of TB following infection (ignoring reinfection) to be 10%
# then calculate the rate of progression from early latency
L_duration <- 0.11
kr_duration <- (((c+u)*((er_duration*L_duration)+(u*L_duration)))-(er_duration*c))/((c+u)*(1-L_duration))

# check it gives us the right lifetime risk and look at how much is in early phase vs late phase
life_risk_duration <- (kr_duration*(c+u) + (er_duration*c))/((er_duration+kr_duration+u)*(c+u))
early_risk_duration <- kr_duration/(kr_duration+er_duration+u)
late_risk_duration <- (er_duration/(kr_duration+er_duration+u))*(c/(c+u))
per_early_duration <- early_risk_duration/life_risk_duration
# calculate the average time spent in the early latent state
time_early_duration <- 1/(er_duration+kr_duration+u)

# Look at cumulative TB assuming no reinfection
t_max <- 100
pars_cum_duration <- c()
Cum_out_duration <- mat.or.vec(21,length(er_duration))
state_cum_duration <- c(LA=1,LB=0,I=0)
for (ii in 1:length(er_duration)){
  
  e <- er_duration[ii]
  k <- kr_duration[ii]
  
  Cum_out_duration[,ii] <- Cum_1(pars_cum_duration,state_cum_duration)[1:21,4]
  
}
Cum_out_duration <- as.data.frame(cbind(seq(0,20),Cum_out_duration))
colnames(Cum_out_duration) <- c("Years",duration_r)
Cum_out_duration <- melt(Cum_out_duration,id.vars="Years")
Cum_plot_duration <- ggplot(data=Cum_out_duration,aes(Years,value))+
  geom_line(aes(colour=variable))+
  ylab("Cumulative proportion with TB")+
  scale_color_discrete(name="Duration of early latency")+
  theme_bw()

# Now look at the transmission model 

# Define beta to run
b_max_duration <- 500
# Array to store required outputs
outputs_duration <- mat.or.vec(b_max_duration*length(er_duration)*length(qr_duration),11)

# counter
kk <- 0
# Pass PT coverage to model - PT is introduced after 600 years (i.e. when model is in equilibrium)
pars_duration <- list(cov_PT=0.05)
# Vector of names for reinfection assumption
reinfection_list <- c("Yes","No")

# time to run model for
tend <- 11

start_time <- Sys.time()
# For each beta 
for (beta in 1:b_max_duration){
  # For each parameterisation
  for (pp in 1:length(er_duration)){
    # For each re-infection assumption
    for (zz in 1:length(qr_duration)){
      
      # counter for run number
      kk <- kk+1
      # Get parameter values
      e <- er_duration[pp]
      k <- kr_duration[pp]
      q <- qr_duration[zz]
      
      # Get steady state and run model 
      source("SS_1.R")
      state <- c(S=S,LA=LA,LB=LB,I=I,PA=0,PB=0,N_PT=0)
      Rag_out <- as.data.frame(TB_model_1(pars_duration,state))
      
      # Store the outputs
      outputs_duration[kk,] <- c(beta,                                                          # beta
                                 100000*Rag_out$Inc[1],                                         # baseline incidence (/100,000)
                                 100*(Rag_out$Inc[1]-Rag_out$Inc[tend])/Rag_out$Inc[1],         # reduction in incidence after 10 years
                                 Rag_out$Number_PT[tend],                                       # PT in year 10
                                 Rag_out$N_PT[tend],                                            # Cumulative PT to year 10
                                 1000000*(Rag_out$Inc[1]-Rag_out$Inc[tend]),                    # Cases averted in year 10 (assuming population of 1,000,000)
                                 1000000*sum(Rag_out$Inc[1]-Rag_out$Inc[1:tend]),               # Cumualive cases averted to year 10 (assuming population of 1,000,000)
                                 Rag_out$Number_PT[tend]/(Rag_out$Inc[tend]-Rag_out$Inc[tend]), # NNT in year 10
                                 Rag_out$N_PT[tend]/sum(Rag_out$Inc[1]-Rag_out$Inc[1:tend]),    # Cumulative NNT to year 10
                                 1/e,                                                           # Duration of early latent
                                 reinfection_list[zz])                                          # Reinfection                                                                       
      
     }
  }
}
end_time <- Sys.time()
end_time - start_time

# Tidy up the outputs           
outputs_duration <- as.data.frame(outputs_duration)               
colnames(outputs_duration) <- c("Beta","Incidence","% Reduction in incidence","Number given PT in year 10","Cumulative number given PT",
                                "Cases averted in year 10 (assuming pop of 1,000,000","Cumulative cases averted",
                                "NNT with PT to avert one case (in year 10)","NNT with PT to avert one case",
                                "Duration of early latency","Reinfection model")

# Plot incidence as a function of beta
outputs_duration_melt <- melt(outputs_duration,id.vars=c("Duration of early latency","Incidence","Reinfection model"))
outputs_duration_melt$value <- as.numeric(as.character(outputs_duration_melt$value))
outputs_duration_melt$Incidence <- as.numeric(as.character(outputs_duration_melt$Incidence))
incidence_plot_duration <- ggplot(outputs_duration_melt[outputs_duration_melt$variable%in%c("Beta"),],aes(value,Incidence))+
  geom_line(aes(colour=as.factor(`Duration of early latency`),linetype=as.factor(`Reinfection model`)),size=1)+
  ylab("TB Incidence /100,000")+
  xlab("Beta")+
  scale_colour_manual(name="Duration of fast latency (yrs)",values=cbPalette)+
  scale_linetype_manual(values=c("solid","dashed"))+
  theme_bw()+theme(legend.position = "bottom")+
  guides(linetype=FALSE)+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

# plot intervention outputs as function of incidence or beta
outputs_duration_melt <- melt(outputs_duration,id.vars=c("Duration of early latency","Incidence","Reinfection model","Beta"))
outputs_duration_melt$value <- as.numeric(as.character(outputs_duration_melt$value))
outputs_duration_melt$Incidence <- as.numeric(as.character(outputs_duration_melt$Incidence))
outputs_duration_melt$Beta <- as.numeric(as.character(outputs_duration_melt$Beta))
outputs_plot_I_duration <- ggplot(outputs_duration_melt[outputs_duration_melt$variable%in%c("% Reduction in incidence",
                                                                                            "NNT with PT to avert one case",
                                                                                            "Cumulative cases averted",
                                                                                            "Cumulative number given PT")&
                                                        outputs_duration_melt$`Reinfection model`=="Yes",],aes(Incidence,value))+
  geom_line(aes(colour=`Duration of early latency`),size=1)+
  facet_wrap(~variable,scales="free",nrow=2)+
  ylab("")+
  xlab("TB Incidence /100,00")+              
  scale_colour_manual(name="Duration of fast latency (yrs)",values=cbPalette)+
  coord_cartesian(xlim = c(-10, 8000))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,NA))+
  scale_x_continuous(expand = c(0, 0))+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(strip.background = element_blank())+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

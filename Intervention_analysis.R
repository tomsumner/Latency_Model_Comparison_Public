# INTERVENTION_ANALYSIS.R

# This script compares the outputs from 3 different model structures using the parameters reported in Menzies and Ragonnet
# It looks at the steady state incidence as a function of the contact parameter
# It simulates a preventive therapy intervention using each model structure and explores different assumptions about the efficacy of therapy

#######################################################################################################################
# SET UP PARAMETERS

# These parameters are the same for all models
tau <- 1         # rate of recovery
u <- 1/50        # natural mortality
m <- 0.03        # TB mortality
qr_S <- c(0.5,0) # rel sus after prior infection - consider 2 situations: 0=no reinfection; 0.5 50% protection against disease following reinfection
tend <- 11       # Time to run to (years)
PT_cov <- 0.05   # annual coverage of PT

# beta values to run to get different baseline incidence - currently set to run values of 1 to 30 in increments of 1 for speed (takes ~ 1 min)
# results presented in paper use values of 0.1 to 30 in increments of 0.1 (seq(0.1,30,by=0.1))
beta_r <- seq(1,30,by=1) 
w_r <- c(0,0.25,0.5)     # values to use for effect of PT (w) - 0=PT provides complete protection against reactivation; 1=PT provides no protection

# These parameters differ by model - need to run "Cumulative_incidence_analysis.R" first to get a_new
# Models 1 and 2 are repeated (in positions 3 and 4) to explore assumptions about number of post PT states (models 1a and 2a)
ar_S <- cbind(c(0,0,0,0,0.0665,a_new[1],0.0665,a_new[1]),c(0,0,0,0,0.085,a_new[2],0.085,a_new[2]))                             
br_S <- cbind(c(0,0.0860,0,0.0860,0,0,0,0),c(0,1-0.91,0,1-0.91,0,0,0,0))                                 
cr_S <- cbind(c(0.000594,0.000594,0.000594,0.000594,0.00337,0.000594,0.00337,0.000594),365*c(5.5e-6,5.5e-6,5.5e-6,5.5e-6,8.52e-6,5.5e-6,8.52e-6,5.5e-6))
kr_S <- cbind(c(0.0826,0.955,0.0826,0.955,0,0,0,0),365*c(1.1e-3,1.1e-2,1.1e-3,1.1e-2,0,0,0,0))      
er_S <- cbind(c(0.872,0,0.872,0,0,0,0,0),365*c(1.1e-2,0,1.1e-2,0,0,0,0,0))                            
                      
#######################################################################################################################
# RUN THE TRANSMISSION MODELS

# Start from steady state and look after 10 years of PT

# Array to store required outputs
outputs_new <- mat.or.vec(length(beta_r)*dim(er_S)[1]*dim(er_S)[2]*length(qr_S)*length(w_r),16)

kk<-0
reinfection_list <- c("Re-infection","No re-infection")
source_list <- c("Menzies","Ragonnet")
model_list <- c("Model 1","Model 2","Model 1","Model 2","Model 3 (original parameters)","Model 3 (updated parameters)","Model 3 (original parameters)","Model 3 (updated parameters)")
PT_model_list <- c("2 states","2 states","1 state","1 state","2 states","2 states","1 state","1 state")
PT_protection_list <- c("100%","75%","50%")

start_time <- Sys.time()

for (pp in 1:dim(er_S)[2]){
  
  for (ii in 1:dim(er_S)[1]){
    
    for (zz in 1:length(qr_S)){
      
      for (bb in 1:length(beta_r)){
        
        for (ww in 1:length(w_r)){
        
          kk <- kk+1
        
          beta <- beta_r[bb] 
        
          a <- ar_S[ii,pp]
          b <- br_S[ii,pp]
          c <- cr_S[ii,pp]
          k <- kr_S[ii,pp]
          e <- er_S[ii,pp]
          q <- qr_S[zz]
          w <- w_r[ww]
        
          # Assume 5% coverage of PT
          pars <- list(cov_PT=PT_cov)
        
          if (ii==1) {
            source("SS_1.R")
            state <- c(S=S,LA=LA,LB=LB,I=I,PA=0,PB=0,N_PT=0)
            Rag_out <- as.data.frame(TB_model_1(pars,state))
          }
          if (ii==2){ 
            source("SS_2.R")
            state <- c(S=S,LA=LA,LB=LB,I=I,PA=0,PB=0,N_PT=0)
            Rag_out <- as.data.frame(TB_model_2(pars,state))
          }
          if (ii==3){ 
            source("SS_1.R")
            state <- c(S=S,LA=LA,LB=LB,I=I,P=0,N_PT=0)
            Rag_out <- as.data.frame(TB_model_1a(pars,state))
          }
          if (ii==4){ 
            source("SS_2.R")
            state <- c(S=S,LA=LA,LB=LB,I=I,P=0,N_PT=0)
            Rag_out <- as.data.frame(TB_model_2a(pars,state))
          }
          if (ii>4) {
            source("SS_3.R")
            state <- c(S=S,L=L,I=I,P=0,N_PT=0)
            Rag_out <- as.data.frame(TB_model_3(pars,state))
          }
        
          outputs_new[kk,] <- c(beta,                                                          # beta
                                model_list[ii],                                                # model
                                source_list[pp],                                               # source
                                reinfection_list[zz],                                          # reinfection  model
                                PT_protection_list[ww],                                        # protection against progression due to PT
                                PT_model_list[ii],                                             # parameterisation
                                Rag_out$Total[tend],                                           # total pop at end to check is 1
                                100000*Rag_out$Inc[1],                                         # baseline incidence (/100,000)
                                100000*Rag_out$Inc[tend],                                      # final incidence (/100,000)
                                100*(Rag_out$Inc[1]-Rag_out$Inc[tend])/Rag_out$Inc[1],         # reduction in incidence after 10 years
                                Rag_out$Number_PT[tend],                                       # PT in year 10
                                10000*Rag_out$N_PT[tend],                                      # Cumulative PT to year 10
                                10000*(Rag_out$Inc[1]-Rag_out$Inc[tend]),                      # Cases averted in year 10 (assuming population of 1,000,000)
                                10000*sum(Rag_out$Inc[1]-Rag_out$Inc[1:tend]),                 # Cumulative cases averted to year 10 (assuming population of 1,000,000)
                                Rag_out$Number_PT[tend]/(Rag_out$Inc[1]-Rag_out$Inc[tend]),    # NNT in year 10
                                Rag_out$N_PT[tend]/sum(Rag_out$Inc[1]-Rag_out$Inc[1:tend]))    # Cumulative NNT to year 10
        
        }
      }
    }
  }
}

end_time <- Sys.time()
end_time - start_time

# Tidy up the outputs
outputs_new<- as.data.frame(outputs_new)               
colnames(outputs_new) <- c("Beta","Model","Source","Reinfection model","PT Protection","PT model","Total","Incidence","Final Incidence",
                         "% Reduction in incidence",
                         "Number given PT in year 10","Cumulative number given PT",
                         "Cases averted in year 10 (assuming pop of 10,000)","Cumulative cases averted",
                         "NNT with PT to avert one case (in year 10)","NNT with PT to avert one case")

outputs_new$Beta = as.numeric(as.character(outputs_new$Beta))
outputs_new$Total = as.numeric(as.character(outputs_new$Total))
outputs_new$Incidence = as.numeric(as.character(outputs_new$Incidence))
outputs_new$`Final Incidence` = as.numeric(as.character(outputs_new$`Final Incidence`))
outputs_new$`% Reduction in incidence` = as.numeric(as.character(outputs_new$`% Reduction in incidence`))
outputs_new$`Number given PT in year 10` = as.numeric(as.character(outputs_new$`Number given PT in year 10`))
outputs_new$`Cumulative number given PT` = as.numeric(as.character(outputs_new$`Cumulative number given PT`))
outputs_new$`Cases averted in year 10 (assuming pop of 10,000)`= as.numeric(as.character(outputs_new$`Cases averted in year 10 (assuming pop of 10,000`))
outputs_new$`Cumulative cases averted` = as.numeric(as.character(outputs_new$`Cumulative cases averted`))
outputs_new$`NNT with PT to avert one case (in year 10)` = as.numeric(as.character(outputs_new$`NNT with PT to avert one case (in year 10)`))
outputs_new$`NNT with PT to avert one case` = as.numeric(as.character(outputs_new$`NNT with PT to avert one case`))

#######################################################################################################################
# GENERATE PLOTS

# Plot steady state incidence as a function of beta - show results for Menzies in the main text
outputs_m2 <- melt(outputs_new,id.vars=c("Beta","Model","Source","Reinfection model","PT model","PT Protection"))

# Figure 3 in main text - steady state incidence for Menzies parameters assuming 2 post PT states
Figure_3 <- ggplot(outputs_m2[which(outputs_m2$variable%in%c("Incidence")&
                         outputs_m2$`PT Protection`=="100%"&
                         outputs_m2$`PT model`=="2 states"&
                         outputs_m2$Source=="Menzies"),],aes(Beta,value))+
                  geom_line(aes(colour=as.factor(Model),linetype=as.factor(`Reinfection model`)),size=1)+
                  ylab("TB Incidence /100,000")+
                  xlab(expression(beta))+
                  scale_colour_manual(name="",values=cbPalette)+
                  scale_linetype_manual(values=c("dashed","solid","dotted"),name="")+
                  scale_y_continuous(expand = c(0, 0))+
                  scale_x_continuous(expand = c(0, 0))+
                  coord_cartesian(ylim = c(0, 1000))+
                  theme_bw()+
                  theme(legend.position="bottom")+
                  theme(strip.background = element_blank())+
                  theme(plot.margin = unit(c(1,1,1,1), "cm"))

# Plot outputs of intervention 
# use which when subsetting as there are NAs in the Incidence column and it returns NA rows unless use the which
outputs_m <- melt(outputs_new[which(outputs_new$Incidence<1200),],id.vars=c("Beta","Model","Source","Reinfection model","Incidence","PT model","PT Protection"))
# Made order of PT protection sensible
outputs_m$PT_f = factor(outputs_m$`PT Protection`, levels=c('0%','25%','50%','75%','100%'))

# Main text figure 4 - Plot Reduction in incidence and NNT as a function of baseline incidence assuming PT protection = 100%
# colours to indicate model, linetype to indicate source

Figure_4 <- ggplot(outputs_m[which(outputs_m$variable%in%c("% Reduction in incidence","NNT with PT to avert one case")&
                                         outputs_m$`Reinfection model`=="Re-infection"&
                                         outputs_m$`PT Protection`=="100%"&
                                         outputs_m$`PT model`=="2 states"),],aes(Incidence,value))+
                geom_line(aes(colour=as.factor(Model),linetype=as.factor(`Source`)),size=1)+
                facet_wrap(~variable,scales="free",nrow=1)+
                xlab("Incidence /100,000")+
                ylab("")+
                scale_colour_manual(name="",values=cbPalette)+
                scale_linetype_manual(values = c("dashed","solid"),name="")+
                coord_cartesian(xlim = c(20, 1000))+
                scale_y_continuous(expand = c(0, 0),limits = c(0,NA))+
                scale_x_continuous(expand = c(0, 0))+
                theme_bw()+theme(legend.position = "bottom")+
                theme(strip.background = element_blank())+
                theme(plot.margin = unit(c(1,1,1,1), "cm"))+
                guides(col = guide_legend(nrow = 1))+
                guides(linetype = guide_legend(nrow = 1))

#######################################################################################################################
# LOOK AT RANGES OF IMPACTS ACROSS MODELS

# Look at range of impacts across "good" model structures (i.e. models 1 and 2)
# Compare to model 3
# See how choice of parameters and assumptions about efficacy of PT affect things

# To do this, we need results at the same incidence - rather than rerunning the model we can interpolate the results we already have

# drop model 3 - plot the results for this as separate lines
x <- outputs_new[outputs_new$Model%in%c("Model 1","Model 2")&outputs_new$`PT Protection`%in%c("100%","75%","50%"),]
x <- droplevels(x)
# reduction in incidence
approx_reduction <- lapply(split(x, list(x$Source,x$Model,x$`Reinfection model`,x$`PT Protection`,x$`PT model`)), function(dat) approx(dat[,8], dat[,10],xout=seq(1,1000,by=1),rule=2)[2])
# NNT
approx_NNT <- lapply(split(x, list(x$Source,x$Model,x$`Reinfection model`,x$`PT Protection`,x$`PT model`)), function(dat) approx(dat[,8], dat[,16],xout=seq(1,1000,by=1),rule=2)[2])

# Now we want to get values at incidence of 1 to 1000 from each approxfun 
interp_all <- c()
for (ii in 1:length(approx_reduction)){
  
  x <- names(approx_reduction[ii])
  temp <- cbind(seq(1,1000,by=1),
                unlist(strsplit(x, ".", fixed = TRUE))[1],
                unlist(strsplit(x, ".", fixed = TRUE))[2],
                unlist(strsplit(x, ".", fixed = TRUE))[3],
                unlist(strsplit(x, ".", fixed = TRUE))[4],
                unlist(strsplit(x, ".", fixed = TRUE))[5],
                unlist(approx_reduction[[ii]]),
                unlist(approx_NNT[[ii]]))
  interp_all <- rbind(interp_all,temp)
  
}

interp_all <- as.data.frame(interp_all)
rownames(interp_all) <- c()
colnames(interp_all) <- c("Incidence","Source","Model","Reinfection","PT","States","Reduction","NNT")

# Calculate ranges in impacts across models for a given source, reinfection model, PT protection and incidence

tt_reduction <- lapply(split(as.numeric(as.character(interp_all$Reduction)),list(interp_all$Source,interp_all$Reinfection,interp_all$PT,interp_all$States,interp_all$Incidence)),summary)
tt_NNT <- lapply(split(as.numeric(as.character(interp_all$NNT)),list(interp_all$Source,interp_all$Reinfection,interp_all$PT,interp_all$States,interp_all$Incidence)),summary)

range_by_source <- c()
for (ii in 1:length(tt_reduction)){
  
  x <- names(tt_reduction[ii])
  temp_reduction <- cbind(unlist(strsplit(x, ".", fixed = TRUE))[1],
                          unlist(strsplit(x, ".", fixed = TRUE))[2],
                          unlist(strsplit(x, ".", fixed = TRUE))[3],
                          unlist(strsplit(x, ".", fixed = TRUE))[4],
                          unlist(strsplit(x, ".", fixed = TRUE))[5],
                          "Reduction",
                          unlist(tt_reduction[[ii]])[1], # min
                          unlist(tt_reduction[[ii]])[4], # mean
                          unlist(tt_reduction[[ii]])[6]) # max
  temp_NNT <- cbind(unlist(strsplit(x, ".", fixed = TRUE))[1],
                    unlist(strsplit(x, ".", fixed = TRUE))[2],
                    unlist(strsplit(x, ".", fixed = TRUE))[3],
                    unlist(strsplit(x, ".", fixed = TRUE))[4],
                    unlist(strsplit(x, ".", fixed = TRUE))[5],
                    "NNT",
                    unlist(tt_NNT[[ii]])[1], # min
                    unlist(tt_NNT[[ii]])[4], # mean
                    unlist(tt_NNT[[ii]])[6]) # max

  range_by_source <- rbind(range_by_source,temp_reduction,temp_NNT)
}
range_by_source <- as.data.frame(range_by_source)
rownames(range_by_source) <- c()
colnames(range_by_source) <- c("Source","Reinfection","PT","States","Incidence","variable","Min","Mean","Max")
    
# Tidy up data - only keep re-infection data
temp<-range_by_source[range_by_source$Reinfection=="Re-infection",]
temp$PT_f = factor(temp$PT, levels=c('100%','75%','50%'))
temp$PT_m = factor(temp$States, levels=c('2 states','1 state'))
temp$Incidence <- as.numeric(as.character(temp$Incidence))
temp$Min <- as.numeric(as.character(temp$Min))
temp$Mean <- as.numeric(as.character(temp$Mean))
temp$Max <- as.numeric(as.character(temp$Max))


# Plot reduction in incidence - figure 5
Figure_5 <- ggplot(temp[temp$variable=="Reduction",],aes(Incidence,Mean))+
                               geom_ribbon(data=temp[temp$variable=="Reduction",],aes(x=Incidence,ymin=Min,ymax=Max,fill=as.factor(Source)),alpha=0.3)+
                               geom_line(data=temp[temp$variable=="Reduction",],aes(Incidence,Min,colour=as.factor(Source)))+
                               geom_line(data=temp[temp$variable=="Reduction",],aes(Incidence,Max,colour=as.factor(Source)))+
                               geom_line(data=outputs_m[outputs_m$variable=="% Reduction in incidence"&
                                              outputs_m$`Reinfection model`=="Re-infection"&
                                              outputs_m$`PT Protection`%in%c("100%","75%","50%")&
                                              outputs_m$Model%in%c("Model 3 (original parameters)","Model 3 (updated parameters)"),],
                               aes(as.numeric(as.character(Incidence)),as.numeric(as.character(value)),colour=as.factor(Source),linetype=as.factor(Model)),size=1)+
                               facet_grid(PT_m~PT_f)+
                               xlab("Incidence /100,000")+
                               ylab("% Reduction in incidence")+
                               scale_colour_manual(name="",values=cbPalette)+
                               scale_fill_manual(name="",values=cbPalette)+
                               scale_linetype_manual(name="",values=c("solid","dashed"))+
                               theme_bw()+theme(legend.position = "bottom")+
                               coord_cartesian(xlim = c(20, 1000))+
                               scale_y_continuous(expand = c(0, 0),limits = c(0,NA))+
                               scale_x_continuous(expand = c(0, 0))+
                               theme(strip.background = element_blank())+
                               theme(plot.margin = unit(c(1,1,1,1), "cm"))+
                               guides(col = guide_legend(nrow = 2))+
                               guides(linetype = guide_legend(nrow = 2))


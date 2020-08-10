# MODEL_3_REPARAMETERISE.R

# This script reparameterises model 3 to give same cumulative incidence as models 1 and 2
# It compares the predictions across models

############################################################################################################
# Calculate new parameters for model 3 to give same lifetime risk as models 1 and 2
# Set c to the values from model 1 and find value of a 
# Add these parameters to the others and run all 4 models

ar_S2 <- rbind(ar_S,(L_S[1,]-(cr_S[1,]/(cr_S[1,]+u)))/(1-(cr_S[1,]/(cr_S[1,]+u))))  # Appendix equation A47
br_S2 <- rbind(br_S,c(0,0))                                 
cr_S2 <- rbind(cr_S,c(0.000594,365*5.5e-6))
kr_S2 <- rbind(kr_S,c(0,0))      
er_S2 <- rbind(er_S,c(0,0))   

###########################################################################################################
# Run the transmission model to calculate the impacts
# Start from steady state and look after 10 years of PT

# Array to store required outputs
outputs_newS <- mat.or.vec(length(beta_r)*dim(er_S2)[1]*dim(er_S2)[2]*length(qr_S)*length(w_r),15)

kk<-0
reinfection_list <- c("Yes","No")
source_list <- c("A","B")
model_list <- c("1","2","3","3 (new parameters)")
PT_protection_list <- c("100% protection","60% protection","25% protection")

start_time <- Sys.time()

for (pp in 1:dim(er_S2)[2]){
  
  for (ii in 1:dim(er_S2)[1]){
    
    for (zz in 1:length(qr_S)){
      
      for (bb in 1:length(beta_r)){
        
        for (ww in 1:length(w_r)){
          
          kk <- kk+1
          
          beta <- beta_r[bb] 
          
          a <- ar_S2[ii,pp]
          b <- br_S2[ii,pp]
          c <- cr_S2[ii,pp]
          k <- kr_S2[ii,pp]
          e <- er_S2[ii,pp]
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
            source("SS_3.R")
            state <- c(S=S,L=L,I=I,P=0,N_PT=0)
            Rag_out <- as.data.frame(TB_model_3(pars,state))
          }
          # Rerun model 3 with reparameterised values
          if (ii==4){   
            source("SS_3.R")
            state <- c(S=S,L=L,I=I,P=0,N_PT=0)
            Rag_out <- as.data.frame(TB_model_3(pars,state))
          }
          
          outputs_newS[kk,] <- c(beta,                                                         # beta
                                model_list[ii],                                                # model
                                source_list[pp],                                               # source
                                reinfection_list[zz],                                          # reinfection  model
                                PT_protection_list[ww],                                        # protection against progression due to PT
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
outputs_newS<- as.data.frame(outputs_newS)               
colnames(outputs_newS) <- c("Beta","Model","Source","Reinfection model","Protection","Total","Incidence","Final Incidence",
                           "% Reduction in TB incidence",
                           "Number given PT in year 10","Cumulative number given PT",
                           "Cases averted in year 10 (assuming pop of 10,000)","Cumulative cases averted",
                           "NNT with PT to avert one case (in year 10)","NNT with PT to avert one case")

outputs_newS$Beta = as.numeric(as.character(outputs_newS$Beta))
outputs_newS$Total = as.numeric(as.character(outputs_newS$Total))
outputs_newS$Incidence = as.numeric(as.character(outputs_newS$Incidence))
outputs_newS$`Final Incidence` = as.numeric(as.character(outputs_newS$`Final Incidence`))
outputs_newS$`% Reduction in TB incidence` = as.numeric(as.character(outputs_newS$`% Reduction in TB incidence`))
outputs_newS$`Number given PT in year 10` = as.numeric(as.character(outputs_newS$`Number given PT in year 10`))
outputs_newS$`Cumulative number given PT` = as.numeric(as.character(outputs_newS$`Cumulative number given PT`))
outputs_newS$`Cases averted in year 10 (assuming pop of 10,000)`= as.numeric(as.character(outputs_newS$`Cases averted in year 10 (assuming pop of 10,000`))
outputs_newS$`Cumulative cases averted` = as.numeric(as.character(outputs_newS$`Cumulative cases averted`))
outputs_newS$`NNT with PT to avert one case (in year 10)` = as.numeric(as.character(outputs_newS$`NNT with PT to avert one case (in year 10)`))
outputs_newS$`NNT with PT to avert one case` = as.numeric(as.character(outputs_newS$`NNT with PT to avert one case`))

#######################################################################################################################
# GENERATE PLOTS
# Plot outputs of intervention including reparameterised version of model 3
# Appendix figure A4
outputs_mS <- melt(outputs_newS[which(outputs_newS$Incidence<1200),],id.vars=c("Beta","Model","Source","Reinfection model","Incidence","Protection"))
outputs_mS$PT_f = factor(outputs_mS$Protection, levels=c('100% protection','60% protection','25% protection'))

Figure_A4 <- ggplot()+
  geom_line(data=outputs_mS[outputs_mS$`Reinfection model`=="Yes"&
                             outputs_mS$variable%in%c('% Reduction in TB incidence'),],
            aes(as.numeric(as.character(Incidence)),as.numeric(as.character(value)),colour=as.factor(Model),linetype=as.factor(Source)),size=1)+
  facet_wrap(~PT_f,scales="fixed",switch="y")+
  xlab("Incidence /100,000")+
  ylab("")+
  scale_colour_manual(name="Model",values=cbPalette)+
  scale_fill_manual(name="Model",values=cbPalette)+
  scale_linetype_manual(values = c("solid","dashed"),name="Parameters")+
  theme_bw()+theme(legend.position = "bottom")+
  coord_cartesian(xlim = c(50, 1000))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,NA))+
  scale_x_continuous(expand = c(0, 0))+
  theme(strip.background = element_blank(),strip.placement = "outside")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  guides(col = guide_legend(nrow = 1))+
  guides(fill = FALSE)+
  guides(linetype = guide_legend(nrow = 1))+
  theme(text = element_text(size=18))


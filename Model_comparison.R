# This code is used to produce the results in "Sumner and White, Modelling the impact of tuberculosis preventive therapy: the importance of model structure"
# It compares the use of different structures for the progression from latency to disease and the impact on the predicted effect of a simple preventive therapy intervention

# Set the working directory
setwd("~/GitHub/Latency_Model_Comparison_Public")

# load libraries
require(deSolve)
require(ggplot2)
require(reshape2)
require(stats)
library(FME)

# load the different models

# these calculate the cumulative TB incidence (ignoring death and re-infection)
source("Model_1_cum.R")
source("Model_2_cum.R")
source("Model_3_cum.R")

# these are the transmission models
source("Model_1.R")
source("Model_2.R")
source("Model_3.R")

# Set up accessible palette for plotting
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##############################################################################################################################################
# This script implements the analysis of the cumulaitve TB risk reported in the appendix
# 
# produces Figure_A1

source("Cumulative_incidence_analysis.R")

##############################################################################################################################################
# This script implements the analysis of the steady state incidence and intervention affects reported in the main text
# 
# produces Figure_2, Figure_3, Figure_A2

source("Intervention_analysis.R")

##############################################################################################################################################
# This script implements the analysis of the duration of the fast latent state for model 1 reported in the appendix
# 
# produces Figure_A3

source("Duration_analysis.R")

##############################################################################################################################################
# This script repeats the intervention analysis with the reparameterised version of model 3
# 
# produces Figure_A4

source("Model_3_reparameterise.R")



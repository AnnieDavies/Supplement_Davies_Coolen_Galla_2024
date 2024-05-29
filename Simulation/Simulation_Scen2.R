# Code to perform the simulation study for scenario 2 in Davies, Coolen and Galla (2024).
# Written by A. Davies (2023).
# Based on code by D. Rizopoulos: https://github.com/drizopoulos/jm_and_lm
# Simulation methodology is described in full in the Supplementary Material of 
# Davies, Coolen and Galla (2024).
# Data is simulated from a joint model with a linear random slopes and random intercept
# longitudinal model and a survival model with a cumulative association.
# Data is split into 10 groups to perform 10-fold cross validation
# The model is fitted & cross validation is performed for the Retarded Kernel models by 
# calling to a Python code using reticulate. The Python results are then transformed to R dataframes
# Then the same analysis is performed for JM and LM in R

library("JMbayes")
library(lcsm)
library("xtable") #for writing to files
library("MASS")
library("splines")
library(tidyverse)
library("reticulate")
library(doSNOW)
library(foreach)

# internal JMbayes function to create landmark data sets
dataLM <- JMbayes:::dataLM

#set wd for reading in functions
setwd("~\\Simulation")
list.files("functions", full.names = TRUE) %>% map(source)

NR <- 50 #number of iterations


## 1. Generate Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(2)
dat_list2 <- list()
for(it in 1:NR){
  print(paste("Data gen iteration ", it))
  dat <- data_gen(scenario = 2, 1000)
  dat_list2 <<- c(dat_list2, list(dat))
}
save.image("DataGenS2.RData")

##2. Split data for cross validation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cl <- makeCluster(detectCores())
registerDoSNOW(cl)

dat_split_list_r <- foreach(it = 1:NR) %dopar% {
  set.seed(it)
  create_split(dat_list2[[it]])
}
parallel::stopCluster(cl)

##3. Delayed Kernel Models~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source_python("DK_parallel.py")

#re-label results so they don't get overwritten
DK.res.scen2 <- DK_results

#export data from python
data_path <- "~\\SAVED DATA\\"
export_DKdata_python(DK.res.scen2, 2, "A", data_path)
export_DKdata_python(DK.res.scen2, 2, "B", data_path)

#convert list to R (to save)
R.DK.results <- convert_py_res(DK.res.scen2)
save.image("DK_S2.RData")

##4. Joint models and Landmarking~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cl <- makeCluster(detectCores(), outfile = "Log-JMLM_Scen2.txt")
registerDoSNOW(cl)

sim.res.scen2 <- foreach(it = 1:NR,
                          .packages = c("JMbayes", "lcsm","xtable", "MASS", "splines")) %dopar% {
                            set.seed(it)
                            run_sim(dat_split_list_r[[it]])
                          }
parallel::stopCluster(cl)

#export results
export_JMdata(sim.res.scen2, scenario = 2, JM = 1, data_path)
export_JMdata(sim.res.scen2, scenario = 2, JM = 2, data_path)
export_LMdata(sim.res.scen2, scenario = 2, data_path)

save.image("JMLM_S2.RData")


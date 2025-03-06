#### Calling program for simulation runs (LRT tests for PBR comparison)
#### Author: Ekkehard Glimm (Aug 2024)
#### QC run by Norbert Hollaender (Sep 2024)  

#Distribute calculations within the loop, then aggregate.

# Load the parallel package
library(clustermq)
source("simulation_functions_myseed.R")
source("analysis_functions_sep_robustvar_fix.R")
source("oneSimRun.R")

# Set the working directory to the directory of the this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Define simulation function
run_simulation_test <- function(scenario, fun_simArm, getPBIRcurve, robustvar, num_cores, num_sims_per_core, alpha) {
  library(haven)
  library(dplyr)
  library(survival)
  scenario_C<-scenario$scenario_C
  scenario_T<-scenario$scenario_T
  scenario_number<-scenario$scenario_number
  myseedlist<-1:num_cores
  num_of_cores<-length(myseedlist)
  num_sims<-num_sims_per_core*num_of_cores
  
  simout<-Q(oneSimRun,myseed=myseedlist, const=list(scenariolist=scenario, reps=num_sims_per_core, 
                                                    fun_simArm=fun_simArm, getPBIRcurve=getPBIRcurve,
                                                    robustvar=robustvar),n_jobs=num_of_cores)
          #oneSimRun does reps repetitions of the sim analysis.
  #aggregation of simulation outputs:
  sim_summary<-data.frame(num_sims=num_sims, power_sand=NA, power_cons=NA,power_lib=NA,
                          d01av=NA,d12av=NA, approxpower01=NA, approxpower12=NA,  
                          Elambda12C=NA, Elambda12T=NA, av_se_robust=NA)
  aggr_list <- do.call(rbind, simout) #aggregate all simulation results into a single long list
  sim_summary$power_sand<-sum(aggr_list$p_sand<alpha)/num_sims
  sim_summary$power_cons<-sum(aggr_list$pval<alpha)/num_sims
  sim_summary$power_lib<-sum(aggr_list$pvalind<alpha)/num_sims
  
  sim_summary$d01av<-sum((aggr_list$numtrans01C+aggr_list$numtrans01T))/num_sims
  sim_summary$approxpower01<-pnorm(qnorm(alpha)-sqrt(sim_summary$d01av/4)*log(scenario_C$lambda01/scenario_T$lambda01))
  sim_summary$d12av<-sum((aggr_list$numtrans12C+aggr_list$numtrans12T))/num_sims
  sim_summary$Elambda12T<-scenario_T$lambda12raw+scenario_T$slope12/(scenario_T$lambda02+scenario_T$lambda01)
  sim_summary$Elambda12C<-scenario_C$lambda12raw+scenario_C$slope12/(scenario_C$lambda02+scenario_C$lambda01)
  sim_summary$approxpower12<-pnorm(qnorm(alpha)-sqrt(sim_summary$d12av/4)*log(sim_summary$Elambda12T/sim_summary$Elambda12C))
  sim_summary$av_se_robust<-sum(sqrt(aggr_list$robustvar))/num_sims
  
  #output:
  simout_all<-list(allruns=aggr_list,summary_sims=sim_summary)
  return(simout_all)
}


#alternative scenario 1:
scenario_C<-data.frame(npatients=300,trtname="C",lambda01=1,lambda02=0.6,lambda12raw=0.5,slope12=0.4,lambda_cens=0.3)
scenario_T<-data.frame(npatients=300,trtname="T",lambda01=1.2,lambda02=0.5,lambda12raw=0.3,slope12=0.2,lambda_cens=0.4)
scenario1 <- list(scenario_C = scenario_C, scenario_T = scenario_T, scenario_number="11_1")

beginsim<-Sys.time()
my1stsimtest<-run_simulation_test(scenario=scenario1, fun_simArm=simulate1trialarm, getPBIRcurve=getPBIRcurve, robustvar=robustvar,
                                  num_cores=10, num_sims_per_core=1000, alpha=0.025)
endsim<-Sys.time()
cat("sim duration with clustermq:", endsim-beginsim, "\n")

#some checks:
#check if all seeds are different:
allmyseeds<-c(my1stsimtest$allruns$seed_C,my1stsimtest$allruns$seed_T)
length(unique(allmyseeds)) #should be 2*num_cores*num_sims_per_core
allmyseeds[duplicated(allmyseeds)]

#write out to file:
# Create a folder for the simulations
if (!dir.exists("SIMRES")) {
  dir.create("SIMRES")
}
# Create a subfolder with the name "scenario_x"
folder_name <- file.path(getwd(), "SIMRES", paste0("ALT_scenario_", scenario1$scenario_number))
if (!dir.exists(folder_name)){dir.create(folder_name)}

# Generate the filenames for the CSV files
file1 <- paste0(folder_name, "/scensetting.csv")
file2 <- paste0(folder_name, "/simoutput.csv")
file3 <- paste0(folder_name, "/sim_summaries.csv")

# Write the data to the CSV files. In all 3 files, first column is the row running index with no column name.
write.csv(scenario1, file1)
write.csv(my1stsimtest$allruns, file2)
write.csv(my1stsimtest$summary_sims, file3)
#### Program: makeplot_Figure2.R 
#### Glimm & Hollaender: Testing PBR (AIMS Mathematics)
#### Generate plots for scenarios A1 to A4 in Figure 2

# Set the working directory to the directory of the this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("analysis_functions_sep_robustvar_fix.R")
source("simulation_functions_myseed.R")

#alternative scenario 1: selected 1st run for Figure 2
num_sims <- 1
scenario_number <- 1
scenario_C<-data.frame(npatients=300,trtname="C",lambda01=1,lambda02=0.6,lambda12raw=0.5,slope12=0.4,lambda_cens=0.3)
scenario_T<-data.frame(npatients=300,trtname="T",lambda01=1.2,lambda02=0.5,lambda12raw=0.3,slope12=0.2,lambda_cens=0.4)

for (rep_sim in 1:num_sims){
  simdata_C<-simulate1trialarm(seed=rep_sim*1239,scenario_C)
  simdata_T<-simulate1trialarm(seed=rep_sim*747,scenario_T)
  simdata<-simdata_C
  simdata$simtimes<-rbind(simdata_C$simtimes,simdata_T$simtimes)
  simdata$eventsdata<-rbind(simdata_C$events,simdata_T$eventsdata)
  
  #PBR curve data:
  output_C<-getPBIRcurve(simdata_C$eventsdata,group='C')
  output_T<-getPBIRcurve(simdata_T$eventsdata,group='T')
  
  title_name <- paste0("A", scenario_number)
  
  plot(stepfun(output_C$stoptime,c(0,output_C$pbin1)),do.points=FALSE,col='red',
       xlim=c(0,12),ylim=c(0,1), main=title_name,
       xlab='time unit',ylab='probability')
  lines(stepfun(output_T$stoptime,c(0,output_T$pbin1)),do.points=FALSE, col='blue')
}  

#alternative scenario 2: selected 7th run for Figure 2
num_sims <- 7
scenario_number <- 2
scenario_C<-data.frame(npatients=300,trtname="C",lambda01=1,lambda02=0.6,lambda12raw=0.5,slope12=0,lambda_cens=0.3)
scenario_T<-data.frame(npatients=300,trtname="T",lambda01=1.5,lambda02=0.6,lambda12raw=0.2,slope12=0.1,lambda_cens=0.3)

for (rep_sim in 1:num_sims){
  simdata_C<-simulate1trialarm(seed=rep_sim*1239,scenario_C)
  simdata_T<-simulate1trialarm(seed=rep_sim*747,scenario_T)
  simdata<-simdata_C
  simdata$simtimes<-rbind(simdata_C$simtimes,simdata_T$simtimes)
  simdata$eventsdata<-rbind(simdata_C$events,simdata_T$eventsdata)
  
  #PBR curve data:
  output_C<-getPBIRcurve(simdata_C$eventsdata,group='C')
  output_T<-getPBIRcurve(simdata_T$eventsdata,group='T')
  
  title_name <- paste0("A", scenario_number)
  plot(stepfun(output_C$stoptime,c(0,output_C$pbin1)),do.points=FALSE,col='red',
       xlim=c(0,12),ylim=c(0,1), main=title_name,
       xlab='time unit',ylab='probability')
  lines(stepfun(output_T$stoptime,c(0,output_T$pbin1)),do.points=FALSE, col='blue')
}  

#alternative scenario 3:  selected 5th run for Figure 2
#num_sims <- 5
#scenario_number <- 3
#scenario_C<-data.frame(npatients=300,trtname="C",lambda01=1.4,lambda02=0.6,lambda12raw=0.7,slope12=0.2,lambda_cens=0.3)
#scenario_T<-data.frame(npatients=300,trtname="T",lambda01=0.9,lambda02=0.6,lambda12raw=0.15,slope12=0.1,lambda_cens=0.3)

# run loop as above


#alternative scenario 4: selected 7th run for Figure 2
#num_sims <- 7
#scenario_number <- 4
#scenario_C<-data.frame(npatients=300,trtname="C",lambda01=0.8,lambda02=0.6,lambda12raw=0.4,slope12=0.3,lambda_cens=0.3)
#scenario_T<-data.frame(npatients=300,trtname="T",lambda01=1.5,lambda02=0.6,lambda12raw=0.7,slope12=0.2,lambda_cens=0.3)

# run loop as above

#Function for simulation run of a simulation scenario on a single core.
#Not intended for direct call, but for being called within calling_program_for_simulation_runs.R.

#Can be called stand-alone, but then seed setting is unnecessarily complicated.
#Calls functions from:
#source("simulation_functions_myseed.R")
#source("analysis_functions_sep_robustvar_fix.R")

#Note on random number seed: The function fun_simArm uses the assigned seed to initiate the random number stream.
#Only within each simulated trial arm, there is a random number stream with the set seed.
#Hence, there are reps streams on every core, each using one seed. 
#To avoid repetition of identical random number streams on different cores, we do this:
#1. Cores get the seed myseed=122+(number of core) 
#2. myseed*rep_sim (control arm) and (myseed+381)*rep_sim (treatment arm) are the seeds used for every instance of
#a simulation within the loop on a core.
#The seeds are all unique:
#2*(myseed*(reps+1)+rep_sim) and 2*(myseed*(reps+1)+rep_sim)-1 obviously produce unique seeds. 
#Added and multiplied different constants to them.
#Probably not necessary, but "masks" pattern in seeds.
oneSimRun<-function(myseed, scenariolist, reps, fun_simArm, getPBIRcurve, robustvar){
  library(haven)
  library(dplyr)
  library(survival)
  source("simulation_functions_myseed.R")
  source("analysis_functions_sep_robustvar_fix.R")
  simoutput<-data.frame(numtrans01C=rep(NA,times=reps), numtrans01T=rep(NA,times=reps),
                        numtrans02C=rep(NA,times=reps), numtrans02T=rep(NA,times=reps),
                        numtrans12C=rep(NA,times=reps), numtrans12T=rep(NA,times=reps),
                        LR01num=rep(NA,times=reps),LR12num=rep(NA,times=reps),
                        LR01denom=rep(NA,times=reps),LR12denom=rep(NA,times=reps),
                        LR=rep(NA,times=reps),pval=rep(NA,times=reps),
                        LRind_num=rep(NA,times=reps),LRind_denom=rep(NA,times=reps),
                        LRind=rep(NA,times=reps),pvalind=rep(NA,times=reps), robustvar=rep(NA,times=reps),
                        LR_sand=rep(NA,times=reps), p_sand=rep(NA,times=reps), seed_C=NA, seed_T=NA)
  for (rep_sim in 1:reps){
    #simulate data:
    simdata_C<-fun_simArm(myseed=(2*(myseed*(reps+1)+rep_sim))*7+381,scenariolist$scenario_C)
    simdata_T<-fun_simArm(myseed=(2*(myseed*(reps+1)+rep_sim)-1)*11+553,scenariolist$scenario_T)
    simdata<-simdata_C
    simdata$simtimes<-rbind(simdata_C$simtimes,simdata_T$simtimes)
    simdata$eventsdata<-rbind(simdata_C$eventsdata,simdata_T$eventsdata)
    
    #PBR curve data:
    output_C<-getPBIRcurve(simdata_C$eventsdata,group='C')
    output_T<-getPBIRcurve(simdata_T$eventsdata,group='T')
    output<-rbind(output_C,output_T)
    output_sorted<-output[order(output$stoptime,output$trt),]
    simoutput$numtrans01C[rep_sim]<-sum((output_sorted$numtrans01)*(output_sorted$trt=="C"))
    simoutput$numtrans01T[rep_sim]<-sum((output_sorted$numtrans01)*(output_sorted$trt=="T"))
    simoutput$numtrans02C[rep_sim]<-sum((output_sorted$numtrans02)*(output_sorted$trt=="C"))
    simoutput$numtrans02T[rep_sim]<-sum((output_sorted$numtrans02)*(output_sorted$trt=="T"))
    simoutput$numtrans12C[rep_sim]<-sum((output_sorted$numtrans12)*(output_sorted$trt=="C"))
    simoutput$numtrans12T[rep_sim]<-sum((output_sorted$numtrans12)*(output_sorted$trt=="T"))
    
    #test:
    simoutput[rep_sim,7:16]<-LRtest(simdata$eventsdata,output_sorted)
    
    #keeping track of distributed computing:
    simoutput$seed_C[rep_sim]<-(2*(myseed*(reps+1)+rep_sim))*7+381
    simoutput$seed_T[rep_sim]<-(2*(myseed*(reps+1)+rep_sim)-1)*11+553
    
    #robust estimates and variance:
    robustvar_out<-robustvar(simdata$eventsdata)
    simoutput$robustvar[rep_sim]<-robustvar_out$robustvar
    #output from robustvar not retained:
      # simout_robustcovar[[rep_sim]]<-robustvar_out$robustcovar
      # simout_indcovar[[rep_sim]]<-robustvar_out$indcovar
      # simout_coxphfit[[rep_sim]]<-robustvar_out$modelout
  } 
  simoutput$LR_sand<-simoutput$LRind_num/sqrt(simoutput$robustvar)
  simoutput$p_sand<-1-pnorm(simoutput$LR_sand)
  return(simoutput)
}
simulate1trialarm<-function(myseed,scenario){
  set.seed(myseed)
  npatients<-scenario$npatients
  trtname<-scenario$trtname
  lambda01<-scenario$lambda01
  lambda02<-scenario$lambda02
  lambda12raw<-scenario$lambda12raw
  slope12<-scenario$slope12
  lambda_cens<-scenario$lambda_cens
  simTimes<-data.frame(trt=rep(trtname,times=npatients),time0x=NA,ev01index=NA,time12=NA,timex2=NA,
                       time_cens=NA, time01obs=NA, time02obs=NA, 
                       totaltime012obs=NA, cens01=NA,cens12=NA)
  #time0x time to transition out of 0
  #ev01index = 1 transition 0->1 otherwise 0->2
  #time12 = time from 1->2
  #cens=Yes or No (Yes=1=TRUE)
  #timex2->time02 if 0->2, time01+time12 if 0->1->2
  simTimes$time0x<-rexp(npatients,rate=lambda01+lambda02)
  simTimes$ev01index<-rbinom(npatients,size=1,prob=lambda01/(lambda01+lambda02))
  
  simTimes$time12<-rexp(npatients,rate=lambda12raw+slope12/lambda01*simTimes$time0x)
  #idea: The duration is response is slowly decreasing(=rate is slowly increasing) as
  #time01 increases
  simTimes$time12<-ifelse(simTimes$ev01index==0,NA,simTimes$time12)
  simTimes$timex2<-ifelse(simTimes$ev01index==0,simTimes$time0x,simTimes$time0x+simTimes$time12)
  
  #independent random censoring:
  simTimes$time_cens<-rexp(npatients,rate=lambda_cens)
  simTimes$cens01<-ifelse(simTimes$time0x<simTimes$time_cens & simTimes$ev01index==1,FALSE,TRUE)
  #FALSE if time 0->1 smaller censoring time
  #TRUE if time 0->1 larger censoring time or if 0->2
  simTimes$cens02<-ifelse(simTimes$time0x<simTimes$time_cens & simTimes$ev01index==0,FALSE,TRUE)
  #FALSE if time 0->2 smaller censoring time
  #TRUE if time 0->2 larger than censoring time or if 0->1
  simTimes$cens12<-ifelse(simTimes$timex2<simTimes$time_cens & simTimes$ev01index==1,FALSE,TRUE)
  #FALSE if total time 0->1->2 smaller censoring time. 
  #TRUE if time 0->1->2 larger than censoring time or if 0->2
  simTimes$time01obs<-ifelse(simTimes$cens01==FALSE, simTimes$time0x,
                             ifelse(simTimes$cens02==TRUE, simTimes$time_cens,simTimes$time0x))
  #time 0->1 if observed
  #time 0->2 if 0->2 observed. In this case cens01=TRUE, so time 0->2 is marked as "censored". I think this is correct, but to be checked.
  #time censored if patient censored in state 0
  simTimes$time02obs<-ifelse(simTimes$cens02==FALSE, simTimes$time0x,
                             ifelse(simTimes$cens01==TRUE, simTimes$time_cens,simTimes$time0x))
  #time 0->2 if observed
  #time 0->1 if 0->1 observed. In this case cens02=TRUE, so time 0->1 is marked as "censored". I think this is correct, but to be checked.
  #time censored if patient censored in state 0
  simTimes$totaltime012obs<-ifelse(simTimes$cens12==FALSE, simTimes$timex2,
                                   ifelse(simTimes$cens01==FALSE, simTimes$time_cens,NA))
  #total time 0->1->2 if observed
  #time censored if patient gets censored in state 1
  #NA if patient gets censored in state 0 or if 0->2 directly 
  
  #Input for analysis is ready:
  #patients in whom 0->1->2 is observed have time01obs, cens01=FALSE, time02obs, cens02=TRUE, totaltime012obs, cens12=FALSE
  #patients in whom 0->1->censored is observed have time01obs, cens01=FALSE, time02obs, cens02=TRUE, totaltime012obs, cens12=TRUE
  #patients in whom 0->2 is observed have time01obs, cens01=TRUE, time02obs, cens02=FALSE, totaltime012obs=NA, cens12=TRUE
  #patients in whom 0->censored is observed have time01obs, cens01=TRUE, time02obs, cens02=TRUE, totaltime012obs=NA, cens12=TRUE
  
  #restructure to create input for test function:
  eventsdata<-data.frame(ID=rep(NA,times=2*nrow(simTimes)),trt=NA,time=NA,trans02=NA,cens0=NA,trans01=NA,trans12=NA,cens1=NA,totcases=NA)
  for (i in 1:nrow(simTimes)){
    if (i==1) count<-i
    #if both times 0->1 and 1->2 are observed:
    if (simTimes$cens01[i]==FALSE & simTimes$cens12[i]==FALSE){
      eventsdata$ID[count]<-i
      eventsdata$trt[count]<-simTimes$trt[i]
      eventsdata$time[count]<-simTimes$time0x[i]
      eventsdata$trans02[count]<-0
      eventsdata$cens0[count]<-0
      eventsdata$trans01[count]<-1
      eventsdata$trans12[count]<-0
      eventsdata$cens1[count]<-0
      eventsdata$totcases[count]<-npatients
      count<-count+1
      eventsdata$ID[count]<-i
      eventsdata$trt[count]<-simTimes$trt[i]
      eventsdata$time[count]<-simTimes$totaltime012obs[i]
      eventsdata$trans02[count]<-0
      eventsdata$cens0[count]<-0
      eventsdata$trans01[count]<-0
      eventsdata$trans12[count]<-1
      eventsdata$cens1[count]<-0
      eventsdata$totcases[count]<-npatients
    }
    #if 0->1 observed, 1->2 censored:
    if (simTimes$cens01[i]==FALSE & simTimes$cens12[i]==TRUE){
      eventsdata$ID[count]<-i
      eventsdata$trt[count]<-simTimes$trt[i]
      eventsdata$time[count]<-simTimes$time0x[i]
      eventsdata$trans02[count]<-0
      eventsdata$cens0[count]<-0
      eventsdata$trans01[count]<-1
      eventsdata$trans12[count]<-0
      eventsdata$cens1[count]<-0
      eventsdata$totcases[count]<-npatients
      count<-count+1
      eventsdata$ID[count]<-i
      eventsdata$trt[count]<-simTimes$trt[i]
      eventsdata$time[count]<-simTimes$totaltime012obs[i]
      eventsdata$trans02[count]<-0
      eventsdata$cens0[count]<-0
      eventsdata$trans01[count]<-0
      eventsdata$trans12[count]<-0
      eventsdata$cens1[count]<-1
      eventsdata$totcases[count]<-npatients
    }
    #if 0->2 observed:
    if (simTimes$cens02[i]==FALSE){
      eventsdata$ID[count]<-i
      eventsdata$trt[count]<-simTimes$trt[i]
      eventsdata$time[count]<-simTimes$time0x[i]
      eventsdata$trans02[count]<-1
      eventsdata$cens0[count]<-0
      eventsdata$trans01[count]<-0
      eventsdata$trans12[count]<-0
      eventsdata$cens1[count]<-0
      eventsdata$totcases[count]<-npatients
    }
    #if censored in 0:
    if (simTimes$cens01[i]==TRUE & simTimes$cens02[i]==TRUE){
      eventsdata$ID[count]<-i
      eventsdata$trt[count]<-simTimes$trt[i]
      eventsdata$time[count]<-simTimes$time_cens[i]
      eventsdata$trans02[count]<-0
      eventsdata$cens0[count]<-1
      eventsdata$trans01[count]<-0
      eventsdata$trans12[count]<-0
      eventsdata$cens1[count]<-0
      eventsdata$totcases[count]<-npatients
    }
    count<-count+1
  }
  eventsdata <- eventsdata[rowSums(is.na(eventsdata)) != ncol(eventsdata), ]
  out<-list(simtimes=simTimes,eventsdata=eventsdata)
  return(out)
}



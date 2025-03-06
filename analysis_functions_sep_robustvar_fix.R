#Functions for the analysis of simulated data in Glimm and Hollaender (2025): Testing probability of being in response.

library(haven)
library(dplyr)
library(survival)

#getPBIRcurve: function to produce the PBR curve with confidence bands. This generates a plot-ready dataset, but no plot.
#See makeplot.R for an example of a call with a plot.
#Structure of input data: 
#colnames(eventsdata)<-c('ID','trt','time','trans02','cens0','trans01','trans12','cens1','totcases')
getPBIRcurve<-function(inputdata,group){
  rawdata_all<-inputdata[order(inputdata$trt, inputdata$time),]
  rawdata<-filter(rawdata_all,trt==group)
  stoptimes_unique<-unique(rawdata$time)
  #Calculate the PBRF curve at event times
  #init:
  p0vec<-c(1,rep(NA,times=length(stoptimes_unique))) #prob to be in state 0 at t0, t1, ...
  p1vec<-c(0,rep(NA,times=length(stoptimes_unique))) #prob to be in state 1 at t0, t1, ...
  p2vec<-c(0,rep(NA,times=length(stoptimes_unique))) #prob to be in state 2 at t0, t1, ...
  #inits for variance estimation:
  Ep0sq<-c(1,rep(NA,times=length(stoptimes_unique))) #aux for variance estimation: expected value of squared quantities
  Ep1sq<-c(0,rep(NA,times=length(stoptimes_unique))) #aux for variance estimation: expected value of squared quantities
  W12<-c(0,rep(NA,times=length(stoptimes_unique)))  #aux for variance estimation: subtracted square of means
  var0<-c(0,rep(NA,times=length(stoptimes_unique))) #variance of the prob of being in state 0 at t0, t1, ...
  var1<-c(0,rep(NA,times=length(stoptimes_unique))) #variance of the prob of being in state 1 at t0, t1, ...
  
  n<-length(unique(rawdata$ID))
  output<-data.frame(stoptime=NA, r0after=NA, r1after=NA, numtrans01=NA, numtrans02=NA, numtrans12=NA, cens0=NA, cens1=NA, pbin0=NA, pbin1=NA, pbin2=NA, var0=NA, var1=NA)

  #Go through all unique times, including censoring times.
  r0<-n #init
  r1<-0 #init
  for (i in 1:length(stoptimes_unique)){
    d01<-sum((rawdata$trans01==1)*(rawdata$time==stoptimes_unique[i]))
      #transitions from state 0 to state 1 at this time
    d02<-sum((rawdata$trans02==1)*(rawdata$time==stoptimes_unique[i]))
    d12<-sum((rawdata$trans12==1)*(rawdata$time==stoptimes_unique[i]))

    c0<-sum((rawdata$cens0==1)*(rawdata$time==stoptimes_unique[i])) #censorings from state 0 at this time
    c1<-sum((rawdata$cens1==1)*(rawdata$time==stoptimes_unique[i])) #censorings from state 0 at this time

    pout0x<-ifelse(r0==0, 0,(d01+d02)/r0)
    pout01<-ifelse(r0==0, 0, d01/r0)
    pout02<-ifelse(r0==0, 0, d02/r0)
    pout12<-ifelse(r1==0, 0, d12/r1)
    p0vec[i+1]=p0vec[i]*(1-pout0x)
    p1vec[i+1]=p0vec[i]*pout01+p1vec[i]*(1-pout12)
    p2vec[i+1]=p0vec[i]*pout02+p1vec[i]*pout12+p2vec[i]

    #variance estimation:
    Eq0xsq<-ifelse(r0==0, 1, pout0x*(1-pout0x)/r0+(1-pout0x)^2)
    Eq12sq<-ifelse(r1==0, (1-pout12)^2, pout12*(1-pout12)/r1+(1-pout12)^2)
    Ep01sq<-ifelse(r0==0, 0, pout01*(1-pout01)/r0+pout01^2)
    Eq12<-1-pout12
    Ep0sq[i+1]<-Ep0sq[i]*Eq0xsq 
    Ep1sq[i+1]<-Ep1sq[i]*Eq12sq+Ep0sq[i]*Ep01sq 
    W12[i+1]<-W12[i]*Eq12^2+p0vec[i]^2*pout01^2
    var1[i+1]<-Ep1sq[i+1]-W12[i+1] #variance of being in state 1 immediately after event time i
    var0[i+1]<-Ep0sq[i+1]-p0vec[i+1]^2 #variance of being in state 0 immediately after event time i

    #update "at risk" (assumption: all censorings happen at event time+epsilon)
    r0<-r0-d01-d02-c0 #in state 0 immediately after this time
    r1<-r1-d12+d01-c1 #in state 1 immediately after this time
    
    output[i,]<-c(stoptimes_unique[i], r0, r1, d01, d02, d12, c0, c1, p0vec[i+1], p1vec[i+1], p2vec[i+1],var0[i+1], var1[i+1])
  }
  output$stdev1<-sqrt(output$var1)
  output$lowlim<-output$pbin1+qnorm(0.025)*output$stdev1
  output$uplim<-output$pbin1+qnorm(0.975)*output$stdev1
  output$trt<-group
  return(output)
}

#getDiffcurve: Function for the difference of two PBR curves. This generates a plot-ready dataset, but no plot.
#See makeplot.R for an example of a call with a plot.
getDiffcurve<-function(inputdata){
  rawdata<-inputdata[order(inputdata$stoptime, inputdata$trt),]
  treatments<-unique(rawdata$trt)
  stoptimes_unique<-unique(rawdata$stoptime)
  #Calculate the PBRF curve at event times
  #init:
  p0Diff<-c(0,rep(NA,times=length(stoptimes_unique))) #prob to be in state 0 at t0, t1, ...
  p1Diff<-c(0,rep(NA,times=length(stoptimes_unique))) #prob to be in state 1 at t0, t1, ...
  p2Diff<-c(0,rep(NA,times=length(stoptimes_unique))) #prob to be in state 2 at t0, t1, ...
  #inits for variance estimation:
  var_p0Diff<-c(0,rep(NA,times=length(stoptimes_unique))) #aux for variance estimation: expected value of squared quantities
  var_p1Diff<-c(0,rep(NA,times=length(stoptimes_unique))) #aux for variance estimation: expected value of squared quantities

  output<-data.frame(stoptime=NA, p0Diff=NA, p1Diff=NA, p2Diff=NA,var_p0Diff=NA,var_p1Diff=NA)
  
  #Go through all unique times, including censoring times.
  #Important note: This must have a beginning at time 0 with all treatments represented at this time.
  #This could be time 0 with no events at all or time 0+epsilon with patients having an immediate state switch.
  #In the latter case, both treatments need to have a non-NA entry in p0, p1, p2.
  
  index_trt1_prev<-1 #init for last valid event time
  index_trt2_prev<-2 #init for last valid event time
  for (i in 1:length(stoptimes_unique)){
    index_trt1<-ifelse(identical(which(rawdata$stoptime==stoptimes_unique[i] & rawdata$trt==treatments[1]),integer(0)),
                       index_trt1_prev,
                       which(rawdata$stoptime==stoptimes_unique[i] & rawdata$trt==treatments[1]))
    index_trt2<-ifelse(identical(which(rawdata$stoptime==stoptimes_unique[i] & rawdata$trt==treatments[2]),integer(0)),
                       index_trt2_prev,
                       which(rawdata$stoptime==stoptimes_unique[i] & rawdata$trt==treatments[2]))
    p0Diff[i]<-rawdata$pbin0[index_trt2]-rawdata$pbin0[index_trt1]
    p1Diff[i]<-rawdata$pbin1[index_trt2]-rawdata$pbin1[index_trt1]
    p2Diff[i]<-rawdata$pbin2[index_trt2]-rawdata$pbin2[index_trt1]
    var_p0Diff[i]<-rawdata$var0[index_trt2]+rawdata$var0[index_trt1]
    var_p1Diff[i]<-rawdata$var1[index_trt2]+rawdata$var1[index_trt1]
    
    output[i,]<-c(stoptimes_unique[i],p0Diff[i],p1Diff[i],p2Diff[i],var_p0Diff[i],var_p1Diff[i])
    index_trt1_prev<-index_trt1
    index_trt2_prev<-index_trt2
  }
  output$stdevdiff1<-sqrt(output$var_p1Diff)
  output$lowlimDiff1<-output$p1Diff+qnorm(0.025)*output$stdevdiff1
  output$uplimDiff1<-output$p1Diff+qnorm(0.975)*output$stdevdiff1
  return(output)
}
#call: diffcurve<-getDiffcurve(output_sorted)

#LRtest: Function for logrank-like tests to compare the two curves:
#cases in eventsdata: Only one of the indicators can be 1. All others must be 0.
#Let 0=no, 1=yes, Cases:
#trans02=1 cens0=0 trans01=0 trans12=0 cens1=0 -> transition from 0 to 2
#trans02=0 cens0=1 trans01=0 trans12=0 cens1=0 -> patient censored in state 0. No state changes for this person
#trans02=0 cens0=0 trans01=1 trans12=0 cens1=0 -> transition from 0 to 1
#trans02=0 cens0=0 trans01=0 trans12=1 cens1=0 -> transition from 1 to 2. Must have been preceded by a transition from 0 to 1 in this patient.
#trans02=0 cens0=0 trans01=0 trans12=0 cens1=1 -> patient censored in state 1. Must have been preceded by a transition from 0 to 1 in this patient.
LRtest<-function(eventsdata,output_sorted){
  eventinx<-as.integer(eventsdata$trans01+2*eventsdata$trans12)
  cdata<-Surv(time=as.numeric(eventsdata$time),event=eventinx, type="mstate")
  #dataset eventsdata converted to Surv object. 
  #Events=0 is a censoring (in either state 0 or state 1 or a death from state 0)
  #Events=1 is a transition from 0 to 1.
  #Events=2  is a transition from 1 to 2.

  stoptimes_unique<-unique(output_sorted$stoptime)
  #at risk in state 0 or 1: initial values
  r0C<-max(eventsdata$totcases*(eventsdata$trt=="C"))
  r0T<-max(eventsdata$totcases*(eventsdata$trt=="T"))
  r1C<-0
  r1T<-0
  #contributions lists to be filled:
  clist01<-data.frame(stoptime=NA,r=NA,rC=NA,rT=NA,d01=NA,d01T=NA,Ed01T=NA,vard01T=NA,tstat_num=NA)
  ncol(clist01)
  clist12<-data.frame(stoptime=NA,r=NA,rC=NA,rT=NA,d12=NA,d12C=NA,Ed12C=NA,vard12C=NA,tstat_num=NA)
  for (i in (1:length(stoptimes_unique))){
    currentdata<-output_sorted[output_sorted$stoptime==stoptimes_unique[i],]
    d01<-sum(currentdata$numtrans01) #total transitions 0-1 at this time
    if (identical(which(currentdata$trt=="T"), integer(0))){d01T<-0}
    else{d01T<-currentdata$numtrans01[which(currentdata$trt=="T")]} #transitions 0-1 in the treatment group. Returns 0 if no line with T exists.
    r0<-r0T+r0C
    if (r0==0){Ed01T<-0}
    else{Ed01T<-d01*r0T/r0} #expected value of d01T under H0
    if (r0>1){vard01T<-d01*r0C*r0T*(r0-d01)/(r0^2*(r0-1))} #variance of d01T (hypergeometric)
    else{vard01T<-0}
    tstat_num<-d01T-Ed01T #numerator of logrank test stat
    clist01[i,]<-c(stoptimes_unique[i],r0,r0C,r0T,d01,d01T,Ed01T,vard01T,tstat_num)
    if (!identical(which(currentdata$trt=="T"), integer(0))){r0T<-currentdata$r0after[which(currentdata$trt=="T")]} #transitions 0-1 in the treatment group. Returns 0 if no line with T exists.
      #otherwise previous roT is retained.
    if (!identical(which(currentdata$trt=="C"), integer(0))){r0C<-currentdata$r0after[which(currentdata$trt=="C")]} #transitions 0-1 in the C group. Returns 0 if no line with C exists.
      #otherwise previous roC is retained.

    #calculations for the 1-2 transition: Note that this uses C as the event group
    d12<-sum(currentdata$numtrans12) #total transitions 1-2 at this time
    if (identical(which(currentdata$trt=="C"), integer(0))){d12C<-0}
    else{d12C<-currentdata$numtrans12[which(currentdata$trt=="C")]} #transitions 1-2 in the C group. Returns 0 if no line with C exists.
    r1<-r1T+r1C
    if(r1==0){Ed12C<-0}
    else{Ed12C<-d12*r1C/r1} #expected value of d12C under H0
    if (r1>1){vard12C<-d12*r1T*r1C*(r1-d12)/(r1^2*(r1-1))} #variance of d12C (hypergeometric)
    else{vard12C<-0} #r1=1 means only one person left at risk, so any event can happen only to this one person,
        #hence variance is 0 (contrib to the numerator is also 0) 
    tstat_num<-d12C-Ed12C #numerator of logrank test stat
    clist12[i,]<-c(stoptimes_unique[i],r1,r1C,r1T,d12,d12C,Ed12C,vard12C,tstat_num)
    if (!identical(which(currentdata$trt=="T"), integer(0))){r1T<-currentdata$r1after[which(currentdata$trt=="T")]} #transitions 1-2 in the treatment group. Returns 0 if no line with T exists.
      #otherwise previous r1T is retained.
    if (!identical(which(currentdata$trt=="C"), integer(0))){r1C<-currentdata$r1after[which(currentdata$trt=="C")]} #transitions 1-2 in the C group. Returns 0 if no line with C exists.
      #otherwise previous r1C is retained.
  }
  #Calculate the test statistic to compare the two curves:
  LR01num<-sum(clist01$tstat_num,na.rm = TRUE)
  LR12num<-sum(clist12$tstat_num,na.rm = TRUE)
  LR01denom<-sqrt(sum(clist01$vard01T,na.rm = TRUE))
  LR12denom<-sqrt(sum(clist12$vard12C,na.rm = TRUE))
  LR<-(LR01num/LR01denom+LR12num/LR12denom)/2 
  LRind_num<-LR01num+LR12num
  LRind_denom<-sqrt(sum(clist01$vard01T,na.rm = TRUE)+sum(clist12$vard12C,na.rm = TRUE))
  LRind<-LRind_num/LRind_denom
 
  pval<-1-pnorm(LR)
  pvalind<-1-pnorm(LRind)
  LRresult<-c(LR01num,LR12num,LR01denom,LR12denom,LR,pval,LRind_num,LRind_denom,LRind,pvalind)
  names(LRresult)<-c("LR01num","LR12num","LR01denom","LR12denom","LR","pval","LRind_num","LRind_denom","LRind","pvalind")
  return(LRresult)
}

#robustvar: function to calculate the robust variance estimate of the LRrob test in Glimm and Hollaender (2025)
robustvar<-function(eventsdata){
  #restructure simdata$eventsdata:
  eventsdata_mod <- data.frame(id = character(),
                               trt = character(),
                               start = numeric(),
                               stop = numeric(),
                               event = numeric(),
                               evtype = character(),
                               stringsAsFactors = FALSE)
  for (i in 1:nrow(eventsdata)){
    current_line<-eventsdata[i,]
    if(current_line$cens0==1){ #censored in state 0
      new_row1 <- data.frame(id = paste(current_line$trt, current_line$ID),
                             trt = current_line$trt,
                             start = 0,
                             stop = current_line$time,
                             event = 0,
                             evtype = "01",
                             stringsAsFactors = FALSE)
      new_row2 <- data.frame(id = paste(current_line$trt, current_line$ID),
                             trt = current_line$trt,
                             start = 0,
                             stop = current_line$time,
                             event = 0,
                             evtype = "02",
                             stringsAsFactors = FALSE)
      new_rows<-rbind(new_row1,new_row2)
    }
    if(current_line$trans01==1){ #0->1 transition observed
      new_row1 <- data.frame(id = paste(current_line$trt, current_line$ID),
                             trt = current_line$trt,
                             start = 0,
                             stop = current_line$time,
                             event = 1,
                             evtype = "01",
                             stringsAsFactors = FALSE)
      new_row2 <- data.frame(id = paste(current_line$trt, current_line$ID),
                             trt = current_line$trt,
                             start = 0,
                             stop = current_line$time,
                             event = 0,
                             evtype = "02",
                             stringsAsFactors = FALSE)
      new_rows<-rbind(new_row1,new_row2)
    }
    if(current_line$trans02==1){ #0->2 transition observed
      new_row1 <- data.frame(id = paste(current_line$trt, current_line$ID),
                             trt = current_line$trt,
                             start = 0,
                             stop = current_line$time,
                             event = 0,
                             evtype = "01",
                             stringsAsFactors = FALSE)
      new_row2 <- data.frame(id = paste(current_line$trt, current_line$ID),
                             trt = current_line$trt,
                             start = 0,
                             stop = current_line$time,
                             event = 1,
                             evtype = "02",
                             stringsAsFactors = FALSE)
      new_rows<-rbind(new_row1,new_row2)
    }
    if(current_line$cens1==1){ #censored in state 1
      start_index <- which(eventsdata_mod$id == paste(current_line$trt, current_line$ID) & eventsdata_mod$evtype == "01")
      # Retrieve the value of "stop" from the identified row
      new_rows <- data.frame(id = paste(current_line$trt, current_line$ID),
                             trt = current_line$trt,
                             start =eventsdata_mod$stop[start_index],
                             stop = current_line$time,
                             event = 0,
                             evtype = "12",
                             stringsAsFactors = FALSE)
    }
    if(current_line$trans12==1){ #1->2 transition observed
      start_index <- which(eventsdata_mod$id == paste(current_line$trt, current_line$ID) & eventsdata_mod$evtype == "01")
      # Retrieve the value of "stop" from the identified row
      new_rows <- data.frame(id = paste(current_line$trt, current_line$ID),
                             trt = current_line$trt,
                             start =eventsdata_mod$stop[start_index],
                             stop = current_line$time,
                             event = 1,
                             evtype = "12",
                             stringsAsFactors = FALSE)
    }
    eventsdata_mod<-rbind(eventsdata_mod,new_rows)
  }
  
  #Model: Let h(x) be the hazard rate for a patient in treatment group x=0 or 1. 
  #This is parametrized as h_j(x)=h_j0*exp(beta_01*x+beta_02*x*(evtype==02)+beta_12*x*(evtype==12))
  #where j=01, 02, 12.
  
  #NAs outside of baseline:
  eventsdata_mod$trt <- factor(eventsdata_mod$trt, levels=c("T","C"))
  eventsdata_mod$evtype <- factor(eventsdata_mod$evtype, levels=c("01","02","12"))
  modelout<-coxph(Surv(time=start, time2=stop, event=event) ~ trt:strata(evtype),
                  cluster = id, data=eventsdata_mod, control = coxph.control(timefix = FALSE))
  
  #estimates of the parameters (NAs removed:)
  coeffests<-na.omit(modelout$coefficients)
  #robust estimate of covariance of the parameters:
  robustcovar<-modelout$var[rowSums(modelout$var^2)>0,colSums(modelout$var^2)>0]
  modeloutind<-coxph(Surv(time=start, time2=stop, event=event) ~ trt:strata(evtype),
                data=eventsdata_mod, control = coxph.control(timefix = FALSE))
  #estimate of covariance of the parameters under independence:
  indcovar<-modeloutind$var[rowSums(modeloutind$var^2)>0,colSums(modeloutind$var^2)>0]
  var_num_LR<-t(c(1,0,1))%*%solve(robustcovar)%*%c(1,0,1)

  output<-list(var_num_LR,robustcovar, indcovar, modelout)
  names(output)<-c("robustvar","robustcovar","indcovar","modelout")
  return(output)
}



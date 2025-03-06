#makeplot.R: example call for rendering a PBR curve and the difference of 2 PBR curves
#Example uses a study dataset from SAS

# Set the working directory to the directory of the this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("analysis_functions_sep_robustvar_fix.R")
#read in data:
eventsdata <- read_sas("eventsdata.sas7bdat") #must be in current working directory
colnames(eventsdata)<-c('ID','trt','time','trans02','cens0','trans01','trans12','cens1','totcases')

#plot two curves with confidence bands
output_BAT<-getPBIRcurve(eventsdata,group='Best available therapy')
output_RUX<-getPBIRcurve(eventsdata,group='Ruxolitinib 10 mg BID')
output_BAT$trt<-'BAT'
output_RUX$trt<-'RUX10mgBID'
output<-rbind(output_BAT,output_RUX)
#output:
# number of rows = number of unique event times (incl censoring times)
#stoptime= event time
#numtransij= number pats going from state i to state j at stoptime
#censj= number pats being censored from state j immediately after stoptime
#pbinj= percentage of patients in state j immediately after stoptime

plot(stepfun(output_RUX$stoptime,c(0,output_RUX$pbin1)),do.points=FALSE,col='red',
     xlim=c(0,max(output$stoptime)+100),ylim=c(-0.1,1.1), main='probability of being in response',
     xlab='days',ylab='probability')
lines(stepfun(output_BAT$stoptime,c(0,output_BAT$pbin1)),do.points=FALSE, col='blue')
lines(stepfun(output_BAT$stoptime,c(0,output_BAT$lowlim)),do.points=FALSE, col='blue',lty='dotted')
lines(stepfun(output_BAT$stoptime,c(0,output_BAT$uplim)),do.points=FALSE, col='blue',lty='dotted')
lines(stepfun(output_RUX$stoptime,c(0,output_RUX$lowlim)),do.points=FALSE, col='red',lty='dotted')
lines(stepfun(output_RUX$stoptime,c(0,output_RUX$uplim)),do.points=FALSE, col='red',lty='dotted')
text(820,0.23,'BAT',cex=0.8)
text(820,0.53,'RUX',cex=0.8)

#plot difference curve:
output_sorted<-output[order(output$stoptime,output$trt),]
diffcurve<-getDiffcurve(output_sorted)

ylimlow<-min(diffcurve$lowlimDiff1)
ylimup<-max(diffcurve$uplimDiff1)
plot(stepfun(diffcurve$stoptime,c(0,diffcurve$p1Diff)),do.points=FALSE,col='red',
     xlim=c(0,max(diffcurve$stoptime)+100),ylim=c(ylimlow-0.05,ylimup+0.05), main='difference of probability of being in response',
     xlab='days',ylab='probability')
lines(stepfun(diffcurve$stoptime,c(0,diffcurve$lowlimDiff1)),do.points=FALSE, col='red',lty='dotted')
lines(stepfun(diffcurve$stoptime,c(0,diffcurve$uplimDiff1)),do.points=FALSE, col='red',lty='dotted')


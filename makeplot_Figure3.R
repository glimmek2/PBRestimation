#makeplot.R: example call for rendering a PBR curve and the difference of 2 PBR curves
#Example uses a study dataset from SAS

# Set the working directory to the directory of the this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("analysis_functions_sep_robustvar_fix.R")
#read in data:
eventsdata <- read.csv("eventsdata.csv") #must be in current working directory
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
     xlim=c(0,max(output$stoptime)+20),ylim=c(-0.1,1.1), main='probability of being in response',
     xlab='days',ylab='probability')
lines(stepfun(output_BAT$stoptime,c(0,output_BAT$pbin1)),do.points=FALSE, col='blue')
lines(stepfun(output_BAT$stoptime,c(0,output_BAT$lowlim)),do.points=FALSE, col='blue',lty='dotted')
lines(stepfun(output_BAT$stoptime,c(0,output_BAT$uplim)),do.points=FALSE, col='blue',lty='dotted')
lines(stepfun(output_RUX$stoptime,c(0,output_RUX$lowlim)),do.points=FALSE, col='red',lty='dotted')
lines(stepfun(output_RUX$stoptime,c(0,output_RUX$uplim)),do.points=FALSE, col='red',lty='dotted')
text(820,0.23,'BAT',cex=0.8)
text(820,0.53,'RUX',cex=0.8)

# Shading the area between lowlim and uplim for BAT
polygon(c(output_BAT$stoptime, rev(output_BAT$stoptime)), c(output_BAT$lowlim, rev(output_BAT$uplim)),
        col = rgb(0, 0, 1, alpha = 0.2), border = NA)

# Shading the area between lowlim and uplim for RUX
polygon(c(output_RUX$stoptime, rev(output_RUX$stoptime)), c(output_RUX$lowlim, rev(output_RUX$uplim)),
        col = rgb(1, 0, 0, alpha = 0.2), border = NA)

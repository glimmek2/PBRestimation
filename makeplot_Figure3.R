# makeplot_Figure3.R 
# Generate Figure 3 in Glimm & Hollaender: testing probability of being in response (AIMS Mathematics)
# Input: REACH-3 data, use R-codes to generate data for PBRcurves

# Set the working directory to the directory of the this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("analysis_functions_sep_robustvar_fix.R")
#read in data:
eventsdata <- read.csv("eventsdata.csv") #must be in current working directory
colnames(eventsdata)<-c('ID','trt','time','trans02','cens0','trans01','trans12','cens1','totcases')


output_BAT<-getPBIRcurve(eventsdata,group='Best available therapy')
output_RUX<-getPBIRcurve(eventsdata,group='Ruxolitinib 10 mg BID')
output_BAT$trt<-'BAT'
output_RUX$trt<-'RUX10mgBID'
output<-rbind(output_BAT,output_RUX)

#####  plot if time in days    

x1 <- output_BAT$stoptime
y1 <- output_BAT$lowlim
y2 <- output_BAT$uplim

v1 <- output_RUX$stoptime
w1 <- output_RUX$lowlim
w2 <- output_RUX$uplim


plot(stepfun(output_RUX$stoptime,c(0,output_RUX$pbin1)),do.points=FALSE,col='blue',
     #xlim=c(0,max(output$stoptime)+100),ylim=c(-0.1,1.1), main='probability of being in response',
     xlim=c(0,800),ylim=c(0,1), main=' ',
     xlab='Days since randomization',ylab='probability of being in response')
polygon(c(x1,rev(x1)),c(y2,rev(y1)),col=adjustcolor("lightpink", alpha.f=0.5), lty = 0)
polygon(c(v1,rev(v1)),c(w2,rev(w1)),col=adjustcolor("lightblue", alpha.f=0.4), lty = 0)
lines(stepfun(output_BAT$stoptime,c(0,output_BAT$pbin1)),do.points=FALSE, col='red')
lines(stepfun(output_RUX$stoptime,c(0,output_RUX$pbin1)),do.points=FALSE, col='blue')
lines(stepfun(output_BAT$stoptime,c(0,output_BAT$lowlim)),do.points=FALSE, col='red',lty='dotted')
lines(stepfun(output_BAT$stoptime,c(0,output_BAT$uplim)),do.points=FALSE, col='red',lty='dotted')
lines(stepfun(output_RUX$stoptime,c(0,output_RUX$lowlim)),do.points=FALSE, col='blue',lty='dotted')
lines(stepfun(output_RUX$stoptime,c(0,output_RUX$uplim)),do.points=FALSE, col='blue',lty='dotted')
text(110,0.29,'BAT',cex=0.8,col='red')
text(110,0.85,'RUX',cex=0.8,col='blue')

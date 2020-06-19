#estimate R0 for Addis Abeba
#R0  estimation using a method given at https://bmcmedinformdecismak.biomedcentral.com/articles/10.1186/1472-6947-12-147#Sec17

library(R0)
#read data
data_csv<-read.csv(file.choose())
View(data_csv)

#subset data for imported case is less than 15%
dat<- data_csv[which(data_csv$Imported.Cases/data_csv$Daily.Cases<=0.15),]
View(dat)
length(data_csv$new.cases) #number of report dates

#Estimate R0
data(data_csv$new.cases)
check.incid(data_csv$new.cases)

mGT<-generation.time("gamma",c(3.7, 0.85)) #c(5.5,2.8))#3-5
#mGT<-generation.time("weibull", c(75, 53.8))#3-5
time_dependent_R0 <- est.R0.TD(data_csv$new.cases, 
                #import = data_csv$Imported.Cases,
                mGT, 
                begin=1, end=26,
                q = c(0.025,0.975), 
                correct = TRUE, nsim = 100000, checked = FALSE)


est.GT(serial.interval=data_csv$new.cases)


#Plot daily R0
plot(time_dependent_R0)
# 95% CI for Ro
time_dependent_R0$conf.int
mean(time_dependent_R0$conf.int$lower)#1.354, 1.322
mean(time_dependent_R0$conf.int$upper)# 3.85, 3.23
#R0 mean, max, min
mean(time_dependent_R0$R)#1.650807, 1.800

max(time_dependent_R0$R) #5.22,6
min(time_dependent_R0$R)#0,0


# R0 weekly
time_dependent_R0.weekly <- smooth.Rt(time_dependent_R0, 7)
time_dependent_R0.weekly
plot(time_dependent_R0.weekly)



##### simulate R0

###data for R0
View(data_R0)
data_R0<-output_SEI_pII_hRD[,c("time","I")]
data_R0$I_ll<-output_SEI_pII_hRD$I
data_R0$I_ul<-output_SEI_pII_hRD$I


#PLOT infectious cases for each RO
plot(data_R0$time,data_R0$I,
     #ylim = c(0,5e+05),
    # xlim = c(0,200),
     type="n",
     xlab='Time',
     ylab = 'Number of cases', main="No intervention",
     panel.first= grid())

lines(data_R0$time,data_R0$I,lwd=3,col='darkgreen')
lines(data_R0$time,data_R0$I_ll,lwd=3,col='blue')
lines(data_R0$time,data_R0$I_ul,lwd=3,col="red")
axis.Date(1, at = seq(times[1], times[length(times)], by="month"),
          labels= seq(times[1], times[length(times)], by="month"),
          format="%d-%m-%Y", las = 2)

legend('topright',c("Upper limit=3.85","Average=2.21",
                    "Lower limit=1.354"),title = "R0",
       lty=c(1,1,1),lwd=c(3,3,3),
       col=c("red","darkgreen",'blue'))



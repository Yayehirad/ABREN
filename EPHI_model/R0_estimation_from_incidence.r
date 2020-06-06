#estimate R0 for Addis Abeba
#R0  estimation using a method given at https://bmcmedinformdecismak.biomedcentral.com/articles/10.1186/1472-6947-12-147#Sec17

library(R0)
#read data
data_csv<-read.csv("Data for Addis Ababa.csv")
View(data_csv)

#subset data for imported case is less than 15%
dat<- data_csv[which(data_csv$Imported.Cases/data_csv$New.Cases<=0.15),]
View(dat)
length(dat$New.Cases) #number of report dates

#Estimate R0
data(dat$New.Cases)
check.incid(dat$New.Cases)

mGT<-generation.time("gamma", c(3.7, 0.85))#3-5
time_dependent_R0 <- est.R0.TD(dat$New.Cases, 
                #import = data_csv$new_imported_case,
                mGT, 
                begin=1, end=24,
                q = c(0.025,0.975), 
                correct = TRUE, nsim = 10000, checked = FALSE)




#Plot daily R0
plot(time_dependent_R0)
# 95% CI for Ro
time_dependent_R0$conf.int
mean(time_dependent_R0$conf.int$lower)#1.354
mean(time_dependent_R0$conf.int$upper)# 3.85
#R0 mean, max, min
mean(time_dependent_R0$R)#2.21

max(time_dependent_R0$R) #5.22
min(time_dependent_R0$R)#0


# R0 weekly
time_dependent_R0.weekly <- smooth.Rt(time_dependent_R0, 7)
time_dependent_R0.weekly
plot(time_dependent_R0.weekly)




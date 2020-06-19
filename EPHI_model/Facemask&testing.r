### Testing preclinical Facemask 70%
View(data_testing)
data_testing<-output_SEI_pII_hRD[,c("time","I")]
data_testing$I_0.1<-output_SEI_pII_hRD$I
data_testing$D_0.1<-output_SEI_pII_hRD$D

data_testing$I_0.15<-output_SEI_pII_hRD$I
data_testing$D_0.15<-output_SEI_pII_hRD$D

data_testing$I_0.2<-output_SEI_pII_hRD$I
data_testing$D_0.2<-output_SEI_pII_hRD$D

data_testing$I_0.25<-output_SEI_pII_hRD$I
data_testing$D_0.25<-output_SEI_pII_hRD$D
#PLOT infectious cases for each I_p isolated
plot(data_testing$time,data_testing$I,
     #ylim = c(0,2.5e+05),
     #xlim = c(0,200),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="Isolation at pre-clinical stage,facemask=70%",
     panel.first= grid())

lines(data_testing$time,data_testing$I,lwd=3,col='red')
lines(data_testing$time,data_testing$I_0.1,lwd=3,col='orange')
lines(data_testing$time,data_testing$I_0.15,lwd=3,col="magenta")
lines(data_testing$time,data_testing$I_0.25,lwd=3,col="blue")

legend('topright',c("0%","10%",
                    "15%","25%"),
       lty=c(1,1,1,1),lwd=c(3,3,3,3),title="Percentage of isolation",
       col=c("red","orange",'magenta',"blue"))




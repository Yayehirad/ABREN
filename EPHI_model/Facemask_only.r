#Facemask only

#data for Facemask comparison
View(data_comparison)
data_comparison<-output_SEI_pII_hRD[,c("time","I","D")]
data_comparison$I_cm_0.25<-output_SEI_pII_hRD$I
data_comparison$D_cm_0.25<-output_SEI_pII_hRD$D
data_comparison$I_cm_0.5<-output_SEI_pII_hRD$I
data_comparison$D_cm_0.5<-output_SEI_pII_hRD$D
data_comparison$I_cm_0.7<-output_SEI_pII_hRD$I
data_comparison$D_cm_0.7<-output_SEI_pII_hRD$D



#PLOT

plot(data_comparison$time,data_comparison$I,
     #ylim = c(0,4e+06),
     #xlim = c(0,200),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="Facemask complience, efficacy=30%", panel.first= grid())

lines(data_comparison$time,data_comparison$I,lwd=3,col='red')
lines(data_comparison$time,data_comparison$I_cm_0.25,lwd=3,col="orange")
lines(data_comparison$time,data_comparison$I_cm_0.5,lwd=3,col='magenta')
lines(data_comparison$time,data_comparison$I_cm_0.7,lwd=3,col='blue')
legend('topright',c("0%","25%","50%",
                    "70%"),title = "Facemask coverage",
       lty=c(1,1,1,1),lwd=c(3,3,3,3),
       col=c("red",'orange','magenta',"blue"))


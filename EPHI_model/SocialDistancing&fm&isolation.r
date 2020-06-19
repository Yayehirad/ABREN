#Facemask=70%, Preclinical isolation = 25% and social distancing

View(data_socialDistancing)
data_socialDistancing<-output_SEI_pII_hRD[,c("time","I")]
data_socialDistancing$I_0.1<-output_SEI_pII_hRD$I
data_socialDistancing$D_0.1<-output_SEI_pII_hRD$D

data_socialDistancing$I_0.2<-output_SEI_pII_hRD$I
data_socialDistancing$D_0.2<-output_SEI_pII_hRD$D

data_socialDistancing$I_0.3<-output_SEI_pII_hRD$I
data_socialDistancing$D_0.3<-output_SEI_pII_hRD$D

#plot
plot(data_socialDistancing$time,data_socialDistancing$I_0.3,
     #ylim = c(0,2.5e+05),
     #xlim = c(0,200),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="Social distancing =30%, isolation=25%,facemask=70%",
     panel.first= grid())

lines(data_socialDistancing$time,data_socialDistancing$I,lwd=3,col='red')
lines(data_socialDistancing$time,data_socialDistancing$I_0.1,lwd=3,col='magenta')
lines(data_socialDistancing$time,data_socialDistancing$I_0.2,lwd=3,col="blue")
lines(data_socialDistancing$time,data_socialDistancing$I_0.3,lwd=3,col="green")

legend('topright',c("0%","10%",
                    "20%","30%"),
       lty=c(1,1,1,1),lwd=c(3,3,3,3),title="SD effectivness",
       col=c("red",'magenta',"blue","green"))





#### time to dates


#SEIRD model for evaluating different estimates of R0

#Step 1 load packages
library(deSolve)

#define function
SEIRD_fun<- function(current_timepoint, state_values, parameters)
{
  
  #creat state variables (local variables)
  S<-state_values[1] # Susceptible 
  E<-state_values[2] # Exposed
  I<-state_values[3] #  Infectious
  
  R<-state_values[4] # Recovered 
  D<-state_values[5]#Died
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      
      dS <- -beta*S*I/N
      
      dE <- beta*S*I/N- sigma*E
      
      dI <- sigma*E - gamma*I- delta*I
      
      dR <- gamma*I
      dD<- delta*I
      
      
      
      #combine results
      results=c(dS,dE,dI,dR,dD)
      list(results)
    }
  )
}




#compartment durations

incubation_period=5 
infectious_period=7 

#force of infection
R0=2.21 #mean 2.21, range 1.354-3.85
effective_contact_rate_per_infectious_per_time= R0/infectious_period #


beta = effective_contact_rate_per_infectious_per_time


# total outflows from compartments
E_outflow = 1/incubation_period

I_outflow_recovery_death = 1/infectious_period

#progression rates from Exposed to Infectious
sigma = E_outflow 
#mortality
proportion_mortality_from_I=0.05
delta= proportion_mortality_from_I*I_u_outflow_recovery_death
#Recovery rate

gamma = I_outflow_recovery_death-delta




#list parameters

parameter_list<- c(beta=beta, 
                   sigma=sigma, 
                   delta=delta,
                   gamma=gamma)
#compartment population 
N=4e+06
E0= N*0.0 #4000
I0=N*0.00025 #1000

R0=N*0 #200
D0=0
S0= N-(E0+I0+R0+D0)

initial_values <- c(S=S0,E=E0,
                    I=I0,
                    R=R0, D=D0)
times=seq(0, 365, by = 1)

output_SEIRD=as.data.frame(lsoda(initial_values, times, SEIRD_fun, parameter_list))


#Check total population N=4e+06
output_SEIRD$N_t<- output_SEIRD$S + output_SEIRD$E + output_SEIRD$I  + output_SEIRD$R +output_SEIRD$D 
View(output_SEIRD)
plot(output_SEIRD$time,output_SEIRD$N_t,type = 'l')
####DATA
#R0_mean
Data_ro<- output_SEIRD[,c("time","I","D")]
#R0_lower LImit
Data_ro$I_r0_min<-output_SEIRD$I 
Data_ro$D_r0_min<-output_SEIRD$D 
#R0_upperlimit
Data_ro$I_r0_max<-output_SEIRD$I 
Data_ro$D_r0_max<-output_SEIRD$D

View(Data_ro)

#PLOT
#PLOT infectious cases for each RO
plot(Data_ro$time,Data_ro$I,
     ylim = c(0,1e+06),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="Infectious cases",
     panel.first= grid())

lines(Data_ro$time,Data_ro$I,lwd=3,col='darkgreen')
lines(Data_ro$time,Data_ro$I_r0_min,lwd=3,col='blue')
lines(Data_ro$time,Data_ro$I_r0_max,lwd=3,col="red")

legend('topright',c("R0=3.85","R0=2.21",
                    "R0=1.354"),
       lty=c(1,1,1),lwd=c(3,3,3),
       col=c("red","darkgreen",'blue'))


#PLOT death for each RO
plot(Data_ro$time,Data_ro$D,
     ylim = c(0,1e+06),
     type="n",
     xlab='Time in days',
     ylab = 'Number of deaths', main="Death",
     panel.first= grid())

lines(Data_ro$time,Data_ro$D,lwd=3,col='darkgreen')
lines(Data_ro$time,Data_ro$D_r0_min,lwd=3,col='blue')
lines(Data_ro$time,Data_ro$D_r0_max,lwd=3,col="red")

legend('topright',c("R0=3.85","R0=2.21",
                    "R0=1.354"),
       lty=c(1,1,1),lwd=c(3,3,3),
       col=c("red","darkgreen",'blue'))












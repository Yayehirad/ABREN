#SEII_hRD model for mask coverage 
#Step 1 load packages
library(deSolve)

#define function
SEIIRD_fun<- function(current_timepoint, state_values, parameters)
{
  
  #creat state variables (local variables)
  S<-state_values[1] # Susceptible 
  E<-state_values[2] # Exposed
  I<-state_values[3] #  Infectious
  I_h<-state_values[4]
  R<-state_values[5] # Recovered 
  D<-state_values[6]#Died
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      
      dS <- -beta_intervention*S*(I+eta_ih*I_h)/N
      
      dE <- beta_intervention*S*(I+ eta_ih*I_h)/N- sigma_i*E - sigma_ih*E
      
      dI <- sigma_i*E - gamma_i*I- delta_i*I - phi_i*I
      
      dI_h <- sigma_ih*E +phi_i*I - gamma_ih*I_h- delta_ih*I_h
      
      dR <- gamma_i*I + gamma_ih*I_h
      
      dD<- delta_i*I + delta_ih*I_h
      
      
      
      #combine results
      results=c(dS,dE,dI,dI_h, dR,dD)
      list(results)
    }
  )
}




#compartment durations

incubation_period=5 
infectious_period=7 
isolation_stay= 10
# total outflows from compartments
E_outflow = 1/incubation_period

I_outflow_recovery_death = 1/infectious_period
I_h_outflow=1/isolation_stay
#force of infection
R0=2.21
effective_contact_rate_per_infectious_per_time= R0/infectious_period #


beta = effective_contact_rate_per_infectious_per_time

#Intervention 
##1 Face Mask
C_M=0.8 #faceMask coverrage c(0,0.25, 0.5, 0.6, 0.7, 0.8)
E_M=0.3#faceMask efficacy
##2 Social distancing 
S_D=1 # relative reduction in contact rate by social distancing

beta_intervention =beta*S_D*(1-C_M*E_M)

## 3 contact tracing and case finding
proportion_E_isolated=0.0 #contact tracing of exposed group
proportion_I_isolated=0.0 #Case finding of Infectious group


relative_infectiousness_of_isolated_compared_to_non_isolated_infectiopu=0.2
eta_ih=relative_infectiousness_of_isolated_compared_to_non_isolated_infectiopu

sigma_ih = proportion_E_isolated*E_outflow 
sigma_i= E_outflow - sigma_ih
#mortality
proportion_mortality_from_I=0.05
proportion_mortality_from_I_h=0.02

delta_i= proportion_mortality_from_I*I_u_outflow_recovery_death
phi_i= proportion_I_isolated*I_u_outflow_recovery_death

delta_ih =proportion_mortality_from_I_h*I_h_outflow
#Recovery rate

gamma_i = I_outflow_recovery_death - delta_i - phi_i
gamma_ih =I_h_outflow - delta_ih




#list parameters

parameter_list<- c(beta_intervention=beta_intervention, 
                   eta_ih=eta_ih,
                   sigma_i=sigma_i, 
                   delta_i=delta_i,
                   gamma_i=gamma_i,
                   
                   sigma_ih=sigma_ih,
                   delta_ih=delta_ih,
                   gamma_ih=gamma_ih)
#compartment population 
N=4e+06
E0= N*0.0 #4000
I0=N*0.00025 #1000
I_h0=N*0
R0=N*0 #200
D0=0
S0= N-(E0+I0+I_h0+R0+D0)

initial_values <- c(S=S0,E=E0,
                    I=I0,I_h=I_h0,
                    R=R0, D=D0)
times=seq(0, 365, by = 1)

output_SEIIRD=as.data.frame(lsoda(initial_values, times, SEIIRD_fun, parameter_list))


#Check total population N=4e+06
output_SEIIRD$N_t<- output_SEIIRD$S + output_SEIIRD$E + output_SEIIRD$I  +
  output_SEIIRD$I_h+ output_SEIIRD$R +output_SEIIRD$D 
#View dataframe
#View(output_SEIIRD)
#plot population and check it is Constant
plot(output_SEIIRD$time,output_SEIIRD$N_t,type = 'l')


# ONLY Mask
#0% mask 
data_mask<- output_SEIIRD[,c("time","I","D")]

#25 % Mask coverage
data_mask$I_.25<- output_SEIIRD$I
data_mask$D_.25<- output_SEIIRD$D
View(data_mask)

#50% mask coverage
data_mask$I_.5<- output_SEIIRD$I
data_mask$D_.5<- output_SEIIRD$D
#60%
data_mask$I_.6<- output_SEIIRD$I
data_mask$D_.6<- output_SEIIRD$D

#70%
data_mask$I_.7<- output_SEIIRD$I
data_mask$D_.7<- output_SEIIRD$D
#80%

data_mask$I_.8<- output_SEIIRD$I
data_mask$D_.8<- output_SEIIRD$D

# save data as csv
write.csv(data_mask,"FaceMask.csv")

#PLOT







#PLOT infectious cases for each face mask coverage
plot(data_mask$time,data_mask$I,
     ylim = c(0,5e+05),
     type="n",
     xlab='Time in days',
     ylab = 'Number of infcetions', main="Infectious cases",
     panel.first= grid())

lines(data_mask$time,data_mask$I,lwd=3,col='red')
lines(data_mask$time,data_mask$I_.25,lwd=3,col='darkorange2')
lines(data_mask$time,data_mask$I_.5,lwd=3,col="pink")
lines(data_mask$time,data_mask$I_.6,lwd=3,col="magenta")
lines(data_mask$time,data_mask$I_.7,lwd=3,col="blue")
lines(data_mask$time,data_mask$I_.8,lwd=3,col="darkgreen")

legend('topright',c("0%","25%","50%","60%","70%","80%"),
       title = "Facemask coverage",
       lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),
       col=c("red","darkorange2",'pink',"magenta","blue","darkgreen"))



#PLOT death for each face mask coverage
plot(data_mask$time,data_mask$D,
     ylim = c(0,3e+05),
     type="n",
     xlab='Time in days',
     ylab = 'Number of infcetions', main="Death",
     panel.first= grid())

lines(data_mask$time,data_mask$D,lwd=3,col='red')
lines(data_mask$time,data_mask$D_.25,lwd=3,col='darkorange2')
lines(data_mask$time,data_mask$D_.5,lwd=3,col="pink")
lines(data_mask$time,data_mask$D_.6,lwd=3,col="magenta")
lines(data_mask$time,data_mask$D_.7,lwd=3,col="blue")
lines(data_mask$time,data_mask$D_.8,lwd=3,col="darkgreen")

legend('topleft',c("0%","25%","50%","60%","70%","80%"),
       title = "Facemask coverage",
       lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),
       col=c("red","darkorange2",'pink',"magenta","blue","darkgreen"))










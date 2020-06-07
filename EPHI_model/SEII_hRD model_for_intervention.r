#SEII_hRD model for Covid for isolation from Exposedby contact tracing
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
##Face Mask
C_M=0.0 #faceMask coverrage
E_M=0.3#faceMask efficacy
##Social distancing 
S_D=1 # relative reduction in contact rate by social distancing

beta_intervention =beta*S_D*(1-C_M*E_M)

## contact tracing and case finding
proportion_E_isolated=0.0 #contact tracing of exposed group
proportion_I_isolated=0.0 #Case finding of Infectious group


relative_infectiousness_of_isolated_compared_to_non_isolated_infectiopu=0.0
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

parameter_list<- c(beta=beta, eta_ih=eta_ih,
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
View(output_SEIIRD)
#plot population and check it is Constant
plot(output_SEIIRD$time,output_SEIIRD$N_t,type = 'l')

#save data 
write.csv(output_SEIIRD,"SEII_hRD_model.csv")

#PLOT

plot(output_SEIIRD$time,output_SEIIRD$R,
     ylim = c(0,4e+06),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="Types of cases", panel.first= grid())

lines(output_SEIIRD$time,output_SEIIRD$S,lwd=2.5,col='blue')
lines(output_SEIIRD$time,output_SEIIRD$E,lwd=2.5,col='darkorchid1')
lines(output_SEIIRD$time,output_SEIIRD$I,lwd=3,col="red")
lines(output_SEIIRD$time,output_SEIIRD$I_h,lwd=3,col="orange")
lines(output_SEIIRD$time,output_SEIIRD$R,lwd=2.5,col = 'green')
lines(output_SEIIRD$time,output_SEIIRD$D,lwd=2.5,col = 'black')
legend('topright',c("Susceptible","Exposed",
                    "Infectious","Isolated infectious",
                    "Recoverd","Dead"),
       lty=c(1,1,1,1,1),lwd=c(3,3,3,3.3),
       col=c("blue",'darkorchid1',"red","orange","green","black"))














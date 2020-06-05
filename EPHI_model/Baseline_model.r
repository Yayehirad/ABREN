#Numeric simulation
#Step 1 load packages
library(deSolve)

#define function
Model_fun<- function(time, state, parameters) {
  
  
  #creat state variables (local variables)
  S<-state[1] # Susceptible 
  E<-state[2] # Exposed
  I_u<-state[3] # non-isolated infectious
  I_h<-state[4] # Isolated/hospitalised infectious
  R<-state[5] # Recovered 
  D<-state[6]# Dead
  
  with(
   as.list(parameters), #variable names within parameters can be used
   {
      #compute derivative
    
    dS <- - beta* S_D*(1-E_M*C_M)*((I_u+ eta_h*I_h)/N)*S
      
    dE <- beta*S_D*(1-E_M*C_M)*((I_u+eta_h*I_h)/N)*S - sigma*E

    dI_u <- (1-f)*sigma*E - p1*phi_u*I_u + p2*delta_u*I_u + (1-(p1+p2))*gamma_u*I_u
     
    dI_h <-  f*sigma*E + p1*phi_u*I_u - (1-q)*gamma_h*I_h + q*delta_h*I_h
     
    dR <- (1-(p1+p2))*gamma_u*I_u + (1-q)*gamma_h*I_h 
      
    dD <- p2*delta_u*I_u + q*delta_h*I_h 
      
      
      #combine results
      results=c(dS,dE,dI_u,dI_h,dR,dD)
      list(results)
    }
  )
  }

#Parameter values
R0=2.6

#compartment durations

incubation_period=5 
infectious_period=7 
hospital_stay=10

detection_delay=5

#force of infection
effective_contact_rate_per_infectious_per_time= R0/infectious_period #
mask_efficacy = 0.3
relative_infectiousness_of_I_h_compared_to_I_u=0.2


beta = effective_contact_rate_per_infectious_per_time
E_M = mask_efficacy
eta_h = relative_infectiousness_of_I_h_compared_to_I_u

# total outflows from compartments
E_outflow = 1/incubation_period

I_u_outflow_recovery_death = 1/infectious_period

I_u_outflow_detection = 1/detection_delay

I_h_outflow = 1/hospital_stay


 
#proportions of Covid related death


proportion_of_I_u_confirmed_and_isolated=0.4
proportion_mortality_from_I_u= 0.02
proportion_mortality_from_I_h= 0.02

p2=proportion_mortality_from_I_u
q=proportion_mortality_from_I_h

#intervention
proportion_of_E_u_progressing_to_I_h=0.0
proportion_population_wearing_faceMask=0.0
reduction_in_contact_rate_by_social_distancing = 1 
proportion_of_I_u_confirmed_and_isolated=0.0

f=proportion_of_E_u_progressing_to_I_h
C_M      = proportion_population_wearing_faceMask
S_D      = reduction_in_contact_rate_by_social_distancing
p1      =proportion_of_I_u_confirmed_and_isolated
phi_u   = p1*I_u_outflow_detection




#progression rates from Exposed to Infectious
sigma = E_outflow 

##Mortality rates
delta_u =  I_u_outflow_recovery_death

delta_h = I_h_outflow

#Recovery rate

gamma_u = I_u_outflow_recovery_death

gamma_h = I_h_outflow


#list parameters

parameters<- c(beta=beta, E_M=E_M, C_M=C_M,
               p1=p1, p2=p2, q=q, f=f,
                eta_h=eta_h,
                sigma=sigma, phi_u=phi_u, 
                   delta_u=delta_u, delta_h=delta_h,
                   gamma_u=gamma_u, gamma_h=gamma_h)
#compartment population 
N=4e+06
E0= N*0.001 #4000
I_u0=N*0.00025 #1000
I_h0=N*0#0.00025 #1000
R0=N*5e-05 #200
D0=N*2.5e-05#100
S0= N-(E0+I_u0+I_h0+R0+D0)

initial_populations <- c(S=S0,E=E0,
                         I_u=I_u0,I_h=I_h0,
                         R=R0,D=D0)
time=seq(0, 365, by = 1)

output_baseline=as.data.frame(lsoda(initial_populations, time, Model_fun, parameters))
View(output_baseline)

#Check total population N=4e+06
output_baseline$N_t<- output_baseline$S + output_baseline$E + output_baseline$I_u + output_baseline$I_h + output_baseline$R + output_baseline$D
View(output_baseline)
plot(output_baseline$time,output_baseline$N_t,type = 'l')

#PLOT
#presentation tuesday @ 6pm 10mins 
plot(output_baseline$time,output_baseline$R,
     #ylim = c(0,4e+06),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="Types of cases", panel.first= grid())

lines(output_baseline$time,output_baseline$S,lwd=2.5,col='chocolate4')
lines(output_baseline$time,output_baseline$E,lwd=2.5,col='darkorchid1')
lines(output_baseline$time,output_baseline$I_u,lwd=3,col="blue")
lines(output_baseline$time,output_baseline$I_h,lwd=3,col = 'red')

lines(output_baseline$time,output_baseline$R,lwd=2.5,col = 'green')
lines(output_baseline$time,output_baseline$D,lwd=2.5,col = 'black')
legend(0,4e+06,c("Susceptible","Exposed",
                 "Non-isolated Infectious",
                 "Isolated Infectious","Recoverd","Dead"),
       lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),
       col=c('chocolate4','darkorchid1',"blue","red","green",'black'))














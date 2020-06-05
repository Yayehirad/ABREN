#Numeric simulation
#load packages
library(deSolve)

#define function
Model_fun<-function(time, state, parameters)
{
  with(as.list(c(state, parameters)), {
  #creat state variables (local variables)
  #S<-state[1] # Susceptible 
  #E_u<-state[2] #Non-quarantined Exposed
  #E_q<-state[3] #Quarantined Exposed
  #I_u<-state[4] # Quarantined infectious
  #I_h<-state[5] # Isolated/hospitalised infectious
  #I_icu<-state[6] # ICU infectious
  #I_a<-state[7] # Asymptomatic infectious
  #R<-state[8] # Recovered 
  #D<-state[9]# Dead
  
  #with(
 #   as.list(parameters), #variable names within parameters can be used
  #  {
      #compute derivative
      
      
 #lambda=(beta*(1-E_M*C_M)*(I_u+eta_a*I_a+eta_h*I_h)/(N-theta_q*(E_q+I_h+I_icu)))
      
     dS = - (beta* S_D* (1-E_M*C_M) * (I_u+eta_a*I_a+eta_h*I_h+eta_icu*I_icu)/N)*S
      
      
      
      dE_u = (beta*S_D*(1-E_M*C_M)*(I_u+eta_a*I_a+eta_h*I_h+eta_icu*I_icu)/N)*S - 
        (sigma_u+alpha_u)*E_u

      
      dE_q = alpha_u*E_u - 
        sigma_q*E_q
     
      
      dI_u = f*sigma_u*E_u - 
        (gamma_u+phi_u+delta_u)*I_u
      
    
      dI_h = sigma_q*E_q+phi_u*I_u+sigma_a*I_a -
        (gamma_h+nu_h+delta_h)*I_h
      
      dI_a = (1-f)*sigma_u*E_u + sigma_q*E_q - 
        (gamma_a+sigma_a+delta_a)*I_a
      
     
      dI_icu = nu_h*I_h - 
        (gamma_icu + delta_icu)*I_icu
      
      
      dR = gamma_u*I_u+ gamma_h*I_h + gamma_a*I_a + gamma_icu*I_icu
      
      dD = delta_u*I_u + delta_h*I_h + delta_a*I_a + delta_icu*I_icu
      
      
      #combine results
      results=c(dS,dE_u, dE_q,dI_u,dI_h,dI_a,dI_icu,dR,dD)
      list(results)
    }
  )
}
#Parameter values
R0=2.5
duration_infectiousness=10
#force of infection
effective_contact_rate_per_infectious_per_time= R0/duration_infectiousness #
mask_efficacy = 0.5
relative_infectiousness_of_I_a_compared_to_I_u=0.5
relative_infectiousness_of_I_h_compared_to_I_u=0.2
relative_infectiousness_of_I_icu_compared_to_I_u=0.2

beta = effective_contact_rate_per_infectious_per_time
E_M = mask_efficacy
eta_a = relative_infectiousness_of_I_a_compared_to_I_u
eta_h = relative_infectiousness_of_I_h_compared_to_I_u
eta_icu = relative_infectiousness_of_I_icu_compared_to_I_u
  
#compartment durations
incubation_period=5 
infectious_period=7 
hospital_stay=10
icu_period=10 



# total outflows from compartments
total_outflow_E_u = 1/incubation_period
total_outflow_E_q = 1/incubation_period
total_outflow_I_u = 1/infectious_period
total_outflow_I_a = 1/infectious_period
total_outflow_I_h = 1/hospital_stay
total_outflow_I_icu = 1/icu_period

 

proportion_of_E_u_progressing_to_I_u=0.4
f=proportion_of_E_u_progressing_to_I_u
#proportions of Covid related death
proportion_mortality_in_I_u= 0.02
proportion_mortality_in_I_a= 0.02
proportion_mortality_in_I_h= 0.02
proportion_mortality_in_I_icu= 0.02
#ICU
proportion_I_h_requiring_ICU=0.1
nu_h = proportion_I_h_requiring_ICU*total_outflow_I_h
#intervention 
proportion_population_wearing_faceMask=0.5
proportion_nonquaratined_Exposed__isolated_by_contact_tracing=0.1
proportion_Asymptomatic_infectious_isolated_by_contact_tracing=0.2
proportion_Nonquarantined_infectious_isolated_by_detection=0.7
reduction_in_contact_rate_by_social_distancing = 1-0.25 

alpha_u  = proportion_nonquaratined_Exposed__isolated_by_contact_tracing * total_outflow_E_u
phi_u    = proportion_Nonquarantined_infectious_isolated_by_detection * total_outflow_I_u
sigma_a  = proportion_Asymptomatic_infectious_isolated_by_contact_tracing * total_outflow_I_a
C_M      = proportion_population_wearing_faceMask
S_D      = reduction_in_contact_rate_by_social_distancing



#progression rates from Exposed to Infectious
sigma_u = total_outflow_E_u - alpha_u
sigma_q = total_outflow_E_q

##Mortality rates
delta_u =  proportion_mortality_in_I_u*total_outflow_I_u
delta_a = proportion_mortality_in_I_a*total_outflow_I_a
delta_h = proportion_mortality_in_I_h*total_outflow_I_h
delta_icu = proportion_mortality_in_I_icu*total_outflow_I_icu




#Recovery rate

gamma_u = total_outflow_I_u-(delta_u + phi_u)

gamma_a = total_outflow_I_a-(delta_a + sigma_a)

gamma_h = total_outflow_I_h - (delta_h + nu_h)

gamma_icu = total_outflow_I_icu - delta_icu 


parameter_list<- c(beta=beta,E_M=E_M,C_M=C_M,
                   eta_h=eta_h,eta_a=eta_a,eta_icu=eta_icu,
                   sigma_u=sigma_u,sigma_a=sigma_a,sigma_q=sigma_q,
                   alpha_u=alpha_u, phi_u=phi_u, f=f,
                   delta_u=delta_u,delta_h=delta_h,delta_a=delta_a,delta_icu=delta_icu,
                   gamma_u=gamma_u,gamma_h=gamma_h,gamma_a=gamma_a,gamma_icu=gamma_icu,
                   nu_h=nu_h)
#compartment population 
N=100000000

E_u0= 0
E_q0=0
I_u0=500
I_a0=500
I_h0=0
I_icu0=0
R0=0
D0=0
S0=N-(E_u0+E_q0+I_u0+I_a0+I_h0+I_icu0+R0+D0)

initial_populations <- c(S=S0,E_u=E_u0,E_q=E_q0,
                         I_u=I_u0,I_a=I_a0,I_h=I_h0,I_icu=I_icu0,
                         R=R0,D=D0)
time=seq(0, 2000, by = 1)

output_baseline=as.data.frame(lsoda(initial_populations, time, Model_fun, parameter_list))
View(output_baseline)














#SEI_pII_hRD model for Covid for isolation from Presympthomatic infectious by contact tracing and Testing
#Step 1 load packages
library(deSolve)


#force of infection

estimate_values <- function(R0, S_D, C_M, E_M){
  effective_contact_rate_per_infectious_per_time= R0/infectious_period
  beta = effective_contact_rate_per_infectious_per_time
  beta_intervention =beta*S_D*(1-C_M*E_M)
  param_values <- list("beta" = beta, "beta_intervention" = beta_intervention)
  return(param_values)
}

#define function
SEI_pII_hRD_fun_first<- function(current_timepoint, state_values, parameters)
{
  
  #creat state variables (local variables)
  S<-state_values[1] # Susceptible 
  E<-state_values[2] # Exposed
  I_p<-state_values[3]
  I<-state_values[4] #  Infectious
  I_h<-state_values[5] #Pre-symptomathic infectious
  R<-state_values[6] # Recovered 
  D<-state_values[7]#Died
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      
      dS <- -return_values_first$beta_intervention*S*(I+eta_ip*I_p+eta_ih*I_h)/N
      
      dE <- return_values_first$beta_intervention*S*(I+eta_ip*I_p+eta_ih*I_h)/N- epsilon*E 
      
      dI_p<- epsilon*E - sigma_ih*I_p - sigma_i*I_p 
      
      dI <- sigma_i*I_p - gamma_i*I- delta_i*I - phi_i*I
      
      dI_h <- sigma_ih*E +phi_i*I - gamma_ih*I_h- delta_ih*I_h
      
      dR <- gamma_i*I + gamma_ih*I_h
      
      dD<- delta_i*I + delta_ih*I_h
      
      
      
      #combine results
      results=c(dS,dE,dI_p,dI,dI_h, dR,dD)
      list(results)
    }
  )
}


SEI_pII_hRD_fun_second<- function(current_timepoint, state_values, parameters)
{
  
  #creat state variables (local variables)
  S<-state_values[1] # Susceptible 
  E<-state_values[2] # Exposed
  I_p<-state_values[3]
  I<-state_values[4] #  Infectious
  I_h<-state_values[5] #Pre-symptomathic infectious
  R<-state_values[6] # Recovered 
  D<-state_values[7]#Died
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      
      dS <- -return_values_second$beta_intervention*S*(I+eta_ip*I_p+eta_ih*I_h)/N
      
      dE <- return_values_second$beta_intervention*S*(I+eta_ip*I_p+eta_ih*I_h)/N- epsilon*E 
      
      dI_p<- epsilon*E - sigma_ih*I_p - sigma_i*I_p 
      
      dI <- sigma_i*I_p - gamma_i*I- delta_i*I - phi_i*I
      
      dI_h <- sigma_ih*E +phi_i*I - gamma_ih*I_h- delta_ih*I_h
      
      dR <- gamma_i*I + gamma_ih*I_h
      
      dD<- delta_i*I + delta_ih*I_h
      
      
      
      #combine results
      results=c(dS,dE,dI_p,dI,dI_h, dR,dD)
      list(results)
    }
  )
}

#compartment durations

incubation_period=3 
presymptomatic_infectious_period=3
infectious_period=7 
isolation_stay= 10
# total outflows from compartments
E_outflow = 1/incubation_period
I_p_outflow =1/presymptomatic_infectious_period
I_outflow = 1/infectious_period
I_h_outflow=1/isolation_stay


## contact tracing and case finding
proportion_I_p_isolated=0.0 #contact tracing and testing of presymptomatic infectious group
proportion_I_isolated=0.0 #Case finding of Infectious group

relative_infectiousness_of_isolated_compared_to_non_isolated_infectious=0.0
relative_infectiousness_of_presympthomatic_compared_to_non_isolated_infectious=0.5
eta_ih = relative_infectiousness_of_isolated_compared_to_non_isolated_infectious
eta_ip = relative_infectiousness_of_presympthomatic_compared_to_non_isolated_infectious

sigma_ih = proportion_I_p_isolated*I_p_outflow 

phi_i= proportion_I_isolated*I_outflow

#progression
epsilon=E_outflow
sigma_i= I_p_outflow - sigma_ih

#mortality
proportion_mortality_from_I=0.05
proportion_mortality_from_I_h=0.02

delta_i= proportion_mortality_from_I*I_outflow

delta_ih =proportion_mortality_from_I_h*I_h_outflow


#Recovery rate

gamma_i = I_outflow - delta_i - phi_i
gamma_ih =I_h_outflow - delta_ih


arrange_paramater <- function(beta){
  # list params
  parameter_list<- c(beta=beta, eta_ih=eta_ih,eta_ip=eta_ip,
                     
                     sigma_i=sigma_i, 
                     delta_i=delta_i,
                     gamma_i=gamma_i,
                     
                     sigma_ih=sigma_ih,
                     delta_ih=delta_ih,
                     gamma_ih=gamma_ih)
  return(parameter_list)
}

R0=2.21
R01=5


#INTERVENTIONS

##Face Mask
C_M=0.0 #faceMask coverrage
E_M=0.3 #faceMask efficacy
##Social distancing 
S_D=1 # relative reduction in contact rate by social distancing

return_values_first <- estimate_values(R0, S_D, E_M, C_M)
return_values_second <- estimate_values(R01, S_D, E_M, C_M)
#compartment population 
N=4e+06
E0= N*0.0 #4000
I_p0=N*0
I0=N*0.00025 #1000
I_h0=N*0
R0=N*0 #200
D0=0
S0= N-(E0+I_p0+I0+I_h0+R0+D0)

initial_values <- c(S=S0,E=E0,
                    I=I0,I_p=I_p0,I_h=I_h0,
                    R=R0, D=D0)
times=seq(0, 365, by = 1)

output_SEI_PII_HRD <- function(sei, param_list){
  output_SEI_pII_hRD=as.data.frame(lsoda(initial_values, times, sei, param_list))
  #Check total population N=4e+06
  output_SEI_pII_hRD$N_t<- output_SEI_pII_hRD$S + output_SEI_pII_hRD$E + output_SEI_pII_hRD$I_p+ 
    output_SEI_pII_hRD$I  +  output_SEI_pII_hRD$I_h+ output_SEI_pII_hRD$R +output_SEI_pII_hRD$D
  
  #View dataframe
  View(output_SEI_pII_hRD)
  #plot population and check it is Constant
  plot(output_SEI_pII_hRD$time,output_SEI_pII_hRD$N_t,type = 'l')
  #save data 
  write.csv(output_SEI_pII_hRD,"SEI_pII_hRD_model.csv")
  
}

output_SEI_PII_HRD(SEI_pII_hRD_fun_first, arrange_paramater(return_values_first$beta))
output_SEI_PII_HRD(SEI_pII_hRD_fun_second, arrange_paramater(return_values_second$beta))


#PLOT

all_plot <- function(value){
  
  
  lines(value$time,value$S,lwd=2.5,col='blue')
  lines(value$time,value$E,lwd=2.5,col='darkorchid1')
  lines(value$time,value$I-p,lwd=2.5,col='darksalmon')
  lines(value$time,value$I,lwd=3,col="red")
  lines(value$time,value$I_h,lwd=3,col="orange")
  lines(value$time,value$R,lwd=2.5,col = 'green')
  lines(value$time,value$D,lwd=2.5,col = 'black')
  legend('topright',c("Susceptible","Exposed","Pre-symptomatic",
                      "Infectious","Isolated infectious",
                      "Recoverd","Dead"),
         lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),
         col=c("blue",'darkorchid1','darksalmon',"red","orange","green","black"))
}


plot(output_SEI_pII_hRD$time,output_SEI_pII_hRD$R,
     ylim = c(0,4e+06),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="R0 = 2.21", panel.first= grid())
all_plot(output_SEI_pII_hRD)

plot(output_SEI_pII_hRD1$time,output_SEI_pII_hRD1$R,
     ylim = c(0,4e+06),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="R0 = 5", panel.first= grid())
all_plot(output_SEI_pII_hRD1)

#SEI_pI_cI_sI_hI_icuRD model for Covid for isolation from Presympthomatic infectious by contact tracing and Testing
#Step 1 load packages
library(deSolve)

#define function
SEI_pI_cI_sI_hI_icuRD<- function(current_timepoint, state_values, parameters)
{
  
  #creat state variables (local variables)
  S<-state_values[1] # Susceptible 
  E<-state_values[2] # Exposed
  I_p<-state_values[3]#preclinical infectious
  I_s<-state_values[4] #Sub-clinical infectious
  I_c<-state_values[5] #  Clinical infectious
  I_h<-state_values[6] # isolated infectious
  I_icu<-state_values[7] #ICU infectious
  R<-state_values[8] # Recovered 
  D<-state_values[9]#Died
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      
      dS <- -beta_intervention*S*(I_c+eta_ip*I_p+eta_is*I_s+eta_ih*I_h+eta_icu*I_icu)/N
      
      dE <- beta_intervention*S*(I_c+eta_ip*I_p+eta_is*I_s+eta_ih*I_h+eta_icu*I_icu)/N- epsilon*E 
      
      dI_p<- epsilon*E - alpha_ip*I_p - sigma_ic*I_p - sigma_is*I_p
      
      dI_s<- sigma_is*I_p - alpha_is*I_s - gamma_is*I_s
      
      dI_c <- sigma_ic*I_p - gamma_ic*I_c - delta_ic*I_c - phi_ic*I_c 
      
      dI_h <- alpha_ip*I_p + phi_ic*I_c +alpha_is*I_s - gamma_ih*I_h - nu_h*I_h
      
      dI_icu<-  nu_h*I_h - gamma_icu*I_icu - delta_icu*I_icu
      
      dR <- gamma_ic*I_c + gamma_ih*I_h + gamma_icu*I_icu +gamma_is*I_s
      
      dD<- delta_ic*I_c + delta_ih*I_h + delta_icu*I_icu
      
      
      
      #combine results
      results=c(dS,dE,dI_p,dI_c,dI_s,dI_h,dI_icu, dR,dD)
      list(results)
    }
  )
}




#compartment durations

incubation_period=3 
preclinical_infectious_period=3
infectious_period=7 
isolation_stay= 10
ICU_stay=20

# total outflows from compartments
E_outflow = 1/incubation_period
I_p_outflow =1/preclinical_infectious_period
I_c_outflow = 1/infectious_period
I_h_outflow=1/isolation_stay
I_s_outflow=1/infectious_period
I_icu_outflow=1/ICU_stay
#force of infection
R0=2.21 #95% CI 1.35-3.3
effective_contact_rate_per_infectious_per_time= R0/infectious_period #


beta = effective_contact_rate_per_infectious_per_time
###
proportion_Subclinical_cases=0.0


#INTERVENTIONS

##Face Mask
C_M=0.0 #faceMask coverrage
E_M=0.3#faceMask efficacy
##Social distancing 
S_D=1 # relative reduction in contact rate by social distancing

beta_intervention =beta*S_D*(1-C_M*E_M)
#


## contact tracing and case finding
proportion_I_p_isolated=0.0 #contact tracing and testing of presymptomatic infectious group
proportion_I_c_isolated=0.0 #Case finding of Infectious group
proportion_I_s_isolated=0.0 #Testing to find Subclinical cases



relative_infectiousness_of_Preclinical_compared_to_non_isolated_Clinical=0.5
relative_infectiousness_of_Subclinical_compared_to_non_isolated_Clinical=0.5
relative_infectiousness_of_Isolated_compared_to_non_isolated_Clinical=0.0
relative_infectiousness_of_ICU_compared_to_non_isolated_Clinical=0.5

eta_ip = relative_infectiousness_of_presympthomatic_compared_to_non_isolated_infectious
eta_is = relative_infectiousness_of_Subclinical_compared_to_non_isolated_Clinical
eta_ih = relative_infectiousness_of_isolated_compared_to_non_isolated_infectious
eta_icu = relative_infectiousness_of_ICU_compared_to_non_isolated_Clinical

alpha_ip = proportion_I_p_isolated*I_p_outflow 
alpha_is = proportion_I_s_isolated*I_s_outflow 

phi_ic= proportion_I_c_isolated*I_c_outflow

#progression
epsilon=E_outflow

sigma_is = proportion_Subclinical_cases*(I_p_outflow-alpha_ip)
sigma_ic= I_p_outflow - sigma_is

# proportion of critical cases, ICU
proportion_I_h_requiring_ICU=0.08

nu_ih=proportion_I_h_requiring_ICU*I_h_outflow


#mortality
proportion_mortality_from_I=0.05

proportion_mortality_from_I_icu=0.5


delta_ic= proportion_mortality_from_I*(I_c_outflow-phi_ic)

delta_icu =proportion_mortality_from_I_icu*I_icu_outflow

#Recovery rate

gamma_ic = I_c_outflow - delta_ic - phi_ic

gamma_ih =I_h_outflow - nu_ih

gamma_is= I_s_outflow - alpha_is

gamma_icu=I_icu_outflow - delta_icu

#list parameters

parameter_list<- c(beta_intervention=beta_intervention, 
                   eta_ih=eta_ih,eta_ip=eta_ip,eta_is=eta_is,
                   
                   sigma_ic=sigma_ic,sigma_is=sigma_is, 
                   alpha_is=alpha_is, alpha_ip=alpha_ip,
                   phi_ic=phi_ic,
                   
                   delta_ic=delta_ic,delta_icu=delta_icu,
                   
                   gamma_ic=gamma_ic, gamma_ih=gamma_ih,gamma_icu=gamma_icu,gamma_is=gamma_is
                   
                        )
#compartment population 
N=4e+06
E0= N*0.0 #4000
I_p0=N*0
I_c0=N*0.00025 #1000
I_s0=N*0
I_h0=N*0
I_icu0=N*0
R0=N*0 #200
D0=N*0
S0= N-(E0+I_p0+I_c0+I_s0+I_h0+I_icu0+R0+D0)

initial_values <- c(S=S0,E=E0,
                    I_c=I_c0,I_p=I_p0,I_s=I_s0,I_h=I_h0,I_icu=I_icu0,
                    R=R0, D=D0)
times=seq(0, 365, by = 1)

output=as.data.frame(lsoda(initial_values, times, SEI_pI_cI_sI_hI_icuRD, parameter_list))
View(output)

#Check total population N=4e+06
output$N_t<- output$S + output$E + output$I_p+ output$I_c + output$I_s+
  output$I_h+ output$I_icu + output$R +output$D 
#View dataframe
View(output)
#plot population and check it is Constant
plot(output$time,output$N_t,type = 'l')

#save data 
write.csv(output_SEI_pII_hRD,"SEI_pII_hRD_model.csv")

#PLOT

plot(output_SEI_pII_hRD$time,output_SEI_pII_hRD$R,
     ylim = c(0,4e+06),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="Types of cases", panel.first= grid())

lines(output_SEI_pII_hRD$time,output_SEI_pII_hRD$S,lwd=2.5,col='blue')
lines(output_SEI_pII_hRD$time,output_SEI_pII_hRD$E,lwd=2.5,col='darkorchid1')
lines(output_SEI_pII_hRD$time,output_SEI_pII_hRD$I-p,lwd=2.5,col='darksalmon')
lines(output_SEI_pII_hRD$time,output_SEI_pII_hRD$I,lwd=3,col="red")
lines(output_SEI_pII_hRD$time,output_SEI_pII_hRD$I_h,lwd=3,col="orange")
lines(output_SEI_pII_hRD$time,output_SEI_pII_hRD$R,lwd=2.5,col = 'green')
lines(output_SEI_pII_hRD$time,output_SEI_pII_hRD$D,lwd=2.5,col = 'black')
legend('topright',c("Susceptible","Exposed","Pre-symptomatic",
                    "Infectious","Isolated infectious",
                    "Recoverd","Dead"),
       lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),
       col=c("blue",'darkorchid1','darksalmon',"red","orange","green","black"))














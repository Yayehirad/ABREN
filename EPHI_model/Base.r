#load packages
library(deSolve)

#define function
Model_fun<-function(current_timepoint, state_values, parameters)
{
  #creat state variables (local variables)
  S_u<-state_values[1] # susceptible 
  S_q<-state_values[2] #early latenyy
  E_u<-state_values[3] #late latency
  E_q<-state_values[4] #infectious
  I_u<-state_values[5]
  I_h<-state_values[6]
  I_cu<-state_values[7]
  I_a<-state_values[8]
  R<-state_values[9]
  D<-state_values[10]
  
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      
      
      dS_u = psi_q*S_q - ((1-p)+(1-q)*p+q*p)*(beta*(1-E_M*C_M)*(I_u+eta_a*I_a+eta_h*I_h)/(N-theta_q*(E_q+I_h+I_cu)))*S_u 
      
      dS_q = (1-p)*(beta*(1-E_M*C_M)*(I_u+eta_a*I_a+eta_h*I_h)/(N-theta_q*(E_q+I_h+I_cu)))*S_u - (theta_j*(beta*(1-E_M*C_M)*(I_u+eta_a*I_a+eta_h*I_h)/(N-theta_q*(E_q+I_h+I_cu)))+psi_q)*S_q
      
      dE_u = (1-q)*p*(beta*(1-E_M*C_M)*(I_u+eta_a*I_a+eta_h*I_h)/(N-theta_q*(E_q+I_h+I_cu)))*S_u - (sigma_u+alpha_u)*E_u
      
      dE_q = q*p*(beta*(1-E_M*C_M)*(I_u+eta_a*I_a+eta_h*I_h)/(N-theta_q*(E_q+I_h+I_cu)))*S_u+alpha_u*E_u+ theta_j*(beta*(1-E_M*C_M)*(I_u+eta_a*I_a+eta_h*I_h)/(N-theta_q*(E_q+I_h+I_cu)))*S_q - sigma_q*E_q
      
      dI_u = f1*sigma_u*E_u - (gamma_u+phi_u+delta_u)*I_u
      
      dI_h = f2*sigma_u*E_u+r*sigma_q*E_q+phi_u*I_u+sigma_a*I_a - (gamma_h+nu_h+delta_h)*I_h
      
      dI_a = (1-f1-f2)*sigma_u*E_u + (1-r)*sigma_q*E_q - (gamma_a+sigma_a+delta_a)*I_a
      
      dI_cu = nu_h*I_h - (gamma_cu + delta_cu)*I_cu
      
      dR = gamma_u*I_u+ gamma_h*I_h + gamma_a*I_a + gamma_cu*I_cu
      
      dD = delta_u*I_u + delta_h*I_h + delta_a*I_a + delta_cu*I_cu
      
      
      #combine results
      results=c(dS_u,dS_q,dE_u, dE_q,dI_u,dI_h,dI_a,dI_cu,dR,dD)
      list(results)
    }
  )
}
#Parameters
#lambda=(beta*(1-E_M*C_M)*(I_u+eta_a*I_a+eta_h*I_h)/(N-theta_q*(E_q+I_h+I_icu)))
beta= 10
p=0.76
q= 0.3
f1=0.4
f2=0.2
psi_q = 1/14
E_M = 0.5
C_M = 0.0567
eta_a = 0.5
eta_h = 0.5
theta_q = 1
theta_j = 0.5
sigma_u = 1/5.1
sigma_q = 1/5.1
sigma_a = 1/4
alpha_u = 0.1168
phi_u = 1/5
r = 0.7

delta_u = 0.015 
delta_h = 0.015 
delta_a = 0.0075
delta_cu = 0.0225
gamma_u = 1/10
gamma_h = 1/8
gamma_a = 1/10
gamma_cu = 1/10
nu_h = 0.083
parameter_list<- c(beta=beta,p=p,q=q, psi_q=psi_q,E_M=E_M,C_M=C_M,
                  eta_h=eta_h,eta_a=eta_a,theta_j=theta_j,
                  sigma_u=sigma_u,sigma_a=sigma_a,sigma_q=sigma_q,
                  alpha_u=alpha_u, phi_u=phi_u, r=r,
                  delta_u=delta_u,delta_h=delta_h,delta_a=delta_a,delta_cu=delta_cu,
                  gamma_u=gamma_u,gamma_h=gamma_h,gamma_a=gamma_a,gamma_cu=gamma_cu,
                  nu_h=nu_h)
N=100000000
S_q0=1000
E_u0= 1000
E_q0=1000
I_u0=100
I_a0=100
I_h0=600
I_cu0=0
R0=160
D0=6
S_u0=N-(S_q0+E_u0+E_q0+I_u0+I_a0+I_h0+I_cu0+R0+D0)

initial_populations <- c(S_u=S_u0,S_q=S_q0,E_u=E_u0,E_q=E_q0,
                         I_u=I_u0,I_a=I_a0,I_h=I_h0,I_cu=I_cu0,
                         R0=R0,D0=D0)
time=seq(0, 365, by = 1)

output_baseline=as.data.frame(lsoda(initial_populations, time, Model_fun, parameter_list))
View(output_baseline)

plot(output_baseline$time,output_baseline$I_h,
     #ylim = c(0,8000000),
     type="n",
     xlab='Time in days',
     ylab = 'Number of cases', main="Types of cases",
     panel.first= grid())

#lines(output$time,output$S,lwd=2.5,col='chocolate4')
#lines(output$time,output$E,lwd=2.5,col='darkorchid1')
lines(output_baseline$time,output_baseline$I_h,lwd=3,col="blue")
lines(output$time,output$I_s,lwd=3,col = 'red')
lines(output$time,output$I_q,lwd=3,col = 'green')
#lines(output$time,output$H,lwd=2.5,col = 'cornflower')
#lines(output$time,output$R,lwd=2.5,col = 'blue')
#lines(output$time,output$D,lwd=2.5,col = 'dark')
legend(0,8e+06,c("Asymptomatic cases","Symptomatic cases",
                 "Quarantine case"),lty=c(1,1,1),lwd=c(3,3,3),
       col=c("blue","red","green"))



#hospital 
output$hospital_capacity<-3000000
plot(output$time,output$H,type = 'n',col='red',lwd=2,
     xlab='Time in days',
     ylab = 'Number of hospitalised cases',
     main="Hospitalisation", panel.first= grid())
lines(output$time,output$H,type = 'l',col='red',lwd=2)
lines(output$time,output$hospital_capacity,type = 'l',col="blue",lwd=2)
lines(output1$time,output1$H,type = 'l',col="green",lwd=2)
legend(0,5e+06,c("No-intervention","Hospital capacity","With intervention"),
       lty=c(1,1,1),lwd=c(2.5,2.5,2.5),
       col=c("red","blue","green"))













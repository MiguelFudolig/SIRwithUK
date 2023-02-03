# This is to try and fit Fudolig and Howard's multi-strain 
# SIR model to UK data with multiple strains.

## Simulating the SIR model

library(tidyverse)
library(ggforce)

#create a function that updates the SIR model
# ini = initial proportions (s,i1,i2,v,r1,r2,d)
# x = list of variables 
#(beta, beta', gamma, gamma', mu,nu)
#where beta and beta' are transmission coefficients, 
# gamma and gamma' are recovery rates
#mu is the mortality rate
#nu is vaccination rate
#tolerance is 1e-9 for 0
mod_sir <- function(ini,
                    N=100000,
                    x,
                    timestart=0,
                    timeend=10, #in weeks
                    dt=1,
                    tol=1e-9){
  sim <- unname(ini)
  
  cases <- matrix(nrow=length(seq(timestart,timeend,by=dt)),
                  ncol=length(ini))
  colnames(cases) <- names(ini)
  cases[1,] <- sim
  
  #probability of transitions from susceptible compartment
  for(t in 2:nrow(cases)){
      p_SI <- (1-exp(-(x$beta*sim[2] + 
                       x$betap*sim[3] + x$nu)))*(x$beta*sim[2])/(x$beta*sim[2] + 
                                                         x$betap*sim[3] + x$nu)
      p_SIn <- (1-exp(-(x$beta*sim[2] + 
                          x$betap*sim[3] + x$nu)))*(x$betap*sim[3])/(x$beta*sim[2] + 
                                                            x$betap*sim[3] + x$nu)
      p_SV <-  (1-exp(-(x$beta*sim[2] + 
                          x$betap*sim[3] + x$nu)))*(x$nu)/(x$beta*sim[2] + x$betap*sim[3] + x$nu)
      
      #probability of transitions from first infected compartment
      p_IR <- (1- exp(-(x$gamma+x$mu)))*(x$gamma/(x$gamma+x$mu))
      p_ID <- (1-exp(-(x$gamma+x$mu)))*(x$mu/(x$gamma+x$mu))
      
      #probability of transition from vaccinated compartment
      p_VIn <- 1-exp(-x$betap*sim[3])
      
      #probability of transition from recovered compartment
      p_RIn <- 1-exp(-x$betap*sim[3])
      
      #probability of transitions from second infected compartment
      p_InRn <- (1- exp(-(x$gammap+x$mu)))*(x$gammap/(x$gammap+x$mu))
      p_InD <- (1-exp(-(x$gammap+x$mu)))*(x$mu/(x$gammap+x$mu))
      
      
      #number of transitions
      #susceptible
      
      S_probs <- c(p_SI, p_SIn, p_SV, 1-sum(p_SI, p_SIn, p_SV))
      n_S <- as.numeric(rmultinom(1,round(N*sim[1]),prob=S_probs))
      
      #first strain
      I_probs <- c(p_IR,p_ID, (1-sum(p_IR,p_ID)))
      n_I <- as.numeric(rmultinom(1,round(N*sim[2]),prob=I_probs))
      
      #new strain
      In_probs <- c(p_InRn, p_InD, 1-sum(p_InRn, p_InD))
      n_In <- as.numeric(rmultinom(1,round(N*sim[3]),prob=In_probs))
      
      #vaccinated
      n_V <- rbinom(1,round(N*sim[4]),prob=p_VIn)
      
      #first recovered
      
      n_R <- rbinom(1,round(N*sim[5]),prob = p_RIn)
      
      #update sim
      change <-(1/N) * c(-sum(n_S[1:3]),
                  n_S[1]-sum(n_I),
                  n_S[2]+n_V+n_R - sum(n_In),
                  n_S[3] - n_V,
                  n_I[1] - n_R,
                  n_In[1],
                  n_I[2]+n_In[2])
      
      sim <- (sim + change)
      sim <- sim*(sim > tol)
      cases[t,] <- sim
  }
  cases  |> 
    as.data.frame() |> 
    mutate(time=seq(timestart,timeend,by=1)) |> 
    mutate_at(vars(time),factor)-> cases
  return(cases)
}

#initial recovered is 1%
initial_recovered <- 0.01
#initial infected based on data
initial_infected <- 0.01
#vaccination rate
vaxrate <- 0.3
#assume 1% of infected is from emergent strain
ini <- c(s=1-vaxrate-initial_infected-initial_recovered,i1=0.99*initial_infected,i2=0.01*initial_infected,v=vaxrate,r1=initial_recovered,r2=0,d=0)
x <- list(beta=2, betap=2, gamma=0.7, gammap=0.7,mu=0.0002, nu=0.0005) #with death and vaccination
# x <- list(beta=2, betap=4.5, gamma=0.7, gammap=0.7,mu=0.000, nu=0.0) #without death and vaccination


test_run <- mod_sir(ini=ini,
                    x=x,
                    timeend = 15)

test_run |>   ggplot(aes(x=time, y=i1*100000)) + geom_point() + geom_point(aes(x=time,y=i2*100000,color="blue"))

test_run |> mutate(infected = i1+i2+0.000001,infprop1=i1/infected,infprop2=i2/infected) |>
  ggplot(aes(x=time)) + #x axis
  geom_point(aes(y=infprop1),size=2) + geom_line(aes(y=infprop1),size=2) + #original
  geom_point(aes(y=infprop2),color="blue", size=2) + geom_line(aes(y=infprop2),color="blue",size=2) +
  geom_hline(yintercept=0.5,linetype="dashed", color="red") + 
  theme_bw() + ylab("Percentage of Infected Population") + xlab("Week after emergence")


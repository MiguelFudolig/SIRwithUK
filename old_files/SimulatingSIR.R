# This is to try and fit Fudolig and Howard's multi-strain 
# SIR model to UK data with multiple strains.

## Simulating the SIR model

library(tidyverse)
library(ggforce)
library(car)
library(emmeans)
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
      sim <- N*sim
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
      
      sim <- (sim/N + change)
      sim <- sim*(sim > tol)
      cases[t,] <- sim
  }
  cases  |> 
    as.data.frame() |> 
    mutate(time=seq(timestart,timeend,by=1)) |> 
    mutate(infected = i1+i2+0.000001,infprop1=i1/infected,infprop2=i2/infected) -> cases 
  
  return(cases)
}


#TESTING ONE RUN=================
#initial recovered is 1%
initial_recovered <- 0.01
#initial infected based on data
initial_infected <- 0.01
#vaccination rate
vaxrate <- 0
#assume 1% of infected is from emergent strain
ini <- c(s=1-vaxrate-initial_infected-initial_recovered,i1=0.99*initial_infected,i2=0.01*initial_infected,v=vaxrate,r1=initial_recovered,r2=0,d=0)
x <- list(beta=1.4, betap=1.25*1.4, gamma=0.7, gammap=0.7,mu=0.0002, nu=0.0005) #with death and vaccination
# x <- list(beta=2, betap=4.5, gamma=0.7, gammap=0.7,mu=0.000, nu=0.0) #without death and vaccination


test_run <- tryCatch({mod_sir(ini=ini,
                              x=x,
                              timestart=0,
                              timeend = 10)},
                     error=function(e){
                       return(data.frame(infprop1=0,infprop2=0))
                     },
                     warning=function(w){
                       return(data.frame(infprop1=0,infprop2=0))
                     }
                     
)

test_run |>   ggplot(aes(x=time, y=i1*100000)) + geom_point() + geom_point(aes(x=time,y=i2*100000,color="blue"))

test_run |> mutate(infected = i1+i2+0.000001,infprop1=i1/infected,infprop2=i2/infected) |>
  ggplot(aes(x=time)) + #x axis
  geom_point(aes(y=infprop1),size=2) + geom_line(aes(y=infprop1),size=2) + #original
  geom_point(aes(y=infprop2),color="blue", size=2) + geom_line(aes(y=infprop2),color="blue",size=2) +
  geom_hline(yintercept=0.5,linetype="dashed", color="red") + 
  theme_bw() + ylab("Percentage of Infected Population") + xlab("Week after emergence") + scale_x_continuous(breaks=0:15)

test_run |> mutate(infected = i1+i2+0.000001,infprop1=log(i1/infected),infprop2=log(i2/infected)) |>
  ggplot(aes(x=time)) + #x axis
  geom_point(aes(y=infprop1),size=2) + geom_line(aes(y=infprop1),size=2) + #original
  geom_point(aes(y=infprop2),color="blue", size=2) + geom_line(aes(y=infprop2),color="blue",size=2) +
  geom_hline(yintercept=0.5,linetype="dashed", color="red") + 
  theme_bw() + ylab("Percentage of Infected Population") + xlab("Week after emergence") + scale_x_continuous(breaks=0:15)

mod1 <- lm(log(infprop2)~time,data=test_run)

coef(mod1)
(log(0.5)- coef(mod1)[1])/coef(mod1)[2]
##LOGISTIC MODELING

# test_run |> mutate(infected = i1+i2+0.000001,infprop1=i1/infected,infprop2=i2/infected) -> test_run

#nonlinear least squares



coef(lm(logit(infprop2)~time,data=test_run)) |> unname()-> gg


logisticmodel <- nls(infprop2~phi1/(1+exp(-(phi2+phi3*time))),
                     start=list(phi1=1,phi2=gg[1], phi3=gg[2]), data=test_run,trace=T)

summary(logisticmodel)

#use coefficients of the logistic model to predict values
#checking if the curve does follow the emergent curve
test_run |> mutate(prediction=predict(logisticmodel,data=test_run)) |> 
  ggplot(aes(x=time)) + #x axis
  geom_point(aes(y=infprop1),size=2) + geom_line(aes(y=infprop1),size=2) + #original
  geom_point(aes(y=infprop2),color="blue", size=2) + geom_line(aes(y=infprop2),color="blue",size=2) +
  geom_line(aes(y=prediction),color="red",size=3)+
  geom_hline(yintercept=0.5,linetype="dashed", color="red") + 
  theme_bw() + ylab("Percentage of Infected Population") + xlab("Week after emergence") + scale_x_continuous(breaks=0:15)


#we use algebra to get the time to 50%.
#t = (-1/phi3) * (phi2 + ln(2*phi1-1))

phi1 <- coef(logisticmodel)[1] |> unname()
phi2 <- coef(logisticmodel)[2] |> unname()
phi3 <- coef(logisticmodel)[3] |> unname()
tdominance <- (-1/phi3) * (phi2 + log(2*phi1-1))

tdominance

#======FUNCTION==============

tdom <- function(ini, x, timestart=0, timeend=10){
  trun <- tryCatch({mod_sir(ini=ini,
                  x=x,
                  timestart=timestart,
                  timeend = timeend)},
                  error=function(e){
                    return(data.frame(infprop1=0,infprop2=0))
                  },
                  warning=function(w){
                    return(data.frame(infprop1=0,infprop2=0))
                  }
                  
  )
  
  if(tail(trun$infprop2,1) <=0.5 | max(trun$infprop2) <=0.5){
    tdominance <- 100 #never dominates
  }
  else{
    #initialize the logistic model using logit transformed model estimates
    coef(lm(logit(infprop2)~time,data=trun)) |> unname()->tcoeffs
    
    phi_2 <- tcoeffs[1]
    phi_3<- tcoeffs[2]
    
    # logisticmodel <- nls(infprop2~phi1/(1+exp(-(phi2+phi3*time))),
    #                      start=list(phi1=1,phi2=phi_2, phi3=phi_3), 
    #                      data=trun,
    #                      trace=F,
    #                      control=list(maxiter=500))
    # 
    # phi1 <- coef(logisticmodel)[1] |> unname()
    # phi2 <- coef(logisticmodel)[2] |> unname()
    # phi3 <- coef(logisticmodel)[3] |> unname()
    # tdominance <- (-1/phi3) * (phi2 + log(2*phi1-1))
    logisticmodel <- tryCatch({nls(infprop2~SSlogis(time,Asym=1,xmid,scal),
                         data=trun,
                         control=list(maxiter=500))},
                         error=function(e){
                           # print(paste("infected probability is ",max(trun$infprop2),", ",tail(trun$infprop2,1)))
                           #rerun with an experimental growth curve
                           rerun <- lm(log(infprop2)~time,
                                        data=trun)
                           
                           return((log(0.5)- coef(rerun)[1])/coef(rerun)[2])
                         },
                         warning=function(w){
                           rerun <- lm(log(infprop2)~time,
                                       data=trun)
                           
                           return((log(0.5)- coef(rerun)[1])/coef(rerun)[2])                         }
                         
                      )
    tdominance <- ifelse(is.numeric(logisticmodel),logisticmodel,coef(logisticmodel)[2])
  }  
  tdominance
}


#testing function
#initial recovered is 1%
initial_recovered <- 0.01
#initial infected based on data
initial_infected <- 0.01
#vaccination rate
vaxrate <- 0.
#assume 1% of infected is from emergent strain
ini <- c(s=1-vaxrate-initial_infected-initial_recovered,i1=0.99*initial_infected,i2=0.01*initial_infected,v=vaxrate,r1=initial_recovered,r2=0,d=0)
x <- list(beta=1.4, betap=1.25*1.75, gamma=0.7, gammap=0.7,mu=0.000, nu=0.000) #with death and vaccination


tdom(ini=ini,x=x,timeend=12)-> tdominance

#=====EXPERIMENT============

#ORIGINAL STRAIN COULD BE gamma =0.7, gammap = 0.7
gamma <- 0.7
gammap <- 0.7
betaset <- c(1.4,1.75,2.1)
betapset <- c(1,1.25,1.5,1.75,2) #multiple of betaset
initial_recovered <- 0.01 #approximation
initial_infected <- 0.01 #loosely based on data
vaxxrateset <- c(0,0.3)
nuset <- c(0, 0.3)
muset <- c(0)
nsims <-100

sim_array <- expand.grid(1:nsims,
                         betaset,
                         betapset,
                         gamma,
                         gammap,
                         vaxxrateset,
                         initial_recovered,
                         initial_infected,
                         nuset,
                         muset) |> as.data.frame()
names(sim_array) <- c("SimNumber",
                      "beta",
                      "betap",
                      "gamma",
                      "gammap",
                      "vaxrate",
                      "ini_recovered",
                      "ini_infected",
                      "nu",
                      "mu")

#experiment function
exp_trun <- function(x1,timestart=0,timeend=12){
  x1 |> unname() |> as.numeric()->x1
  #initial recovered is 1%
  initial_recovered <- x1[7]
  #initial infected based on population
  initial_infected <- x1[8]
  #vaccination rate
  vaxrate <- x1[6]
  
  x <-list(beta=x1[2], 
           betap=x1[3]*x1[2],
           gamma=x1[4], 
           gammap=x1[5],
           nu=x1[9], 
           mu=x1[10]) #no birth and vaccination

  ini <- c(s=1-vaxrate-initial_infected-initial_recovered,
           i1=0.99*initial_infected,
           i2=0.01*initial_infected,
           v=vaxrate,
           r1=initial_recovered,
           r2=0,
           d=0)
  tdom(ini=ini,x=x,timestart=timestart,timeend=timeend)-> tdominance
  tdominance
}



rm(expresults)
sim_array |> rowwise() |> 
  mutate(tdom = exp_trun(c(SimNumber,beta,betap,gamma,gammap,vaxrate,
                           ini_recovered,ini_infected,nu,mu))) |>  ungroup()|> 
  mutate_at(vars(beta,betap,gamma,gammap,vaxrate,
                 ini_recovered,ini_infected,nu,mu),factor)-> expresults


expresults |>  group_by(beta,betap,gamma,gammap,vaxrate,
                               ini_recovered,ini_infected,nu,mu) |>
  mutate(nondom=ifelse(tdom==100,1,0)) |> 
  summarise(mean=mean(tdom),sd=sd(tdom),max=max(tdom),min=min(tdom),
            prop_nondom=sum(nondom)/n(),median=quantile(tdom,0.50)) ->sumstats

sumstats |> filter(mean > 12) |> View()
mod1 <- lm(tdom~beta*betap*vaxrate*nu,data=expresults,contrasts = list(beta=contr.SAS,
                                                                          betap=contr.SAS,
                                                                          nu=contr.SAS))


# mod1 <- lm(tdom~beta*betap*vaxrate*nu*mu,data=expresults,contrasts = list(beta=contr.SAS,
#                                                                betap=contr.SAS,
#                                                                nu=contr.SAS))
mod1 |> aov() |> summary()
emmeans(mod1,~nu*beta*betap*vaxrate)
emmeans(mod1,pairwise~nu|beta*betap*vaxrate,adjust="tukey")

# emmeans(mod1,pairwise~mu|beta*betap*vaxrate*nu,adjust="tukey")


#analyzing without non-dominating runs:

expresults |> filter(tdom !=100) ->expresults_dom



expresults_dom |>  group_by(beta,betap,gamma,gammap,vaxrate,
                        ini_recovered,ini_infected,nu,mu) |>
  summarise(mean=mean(tdom),sd=sd(tdom),max=max(tdom),min=min(tdom),
            median=quantile(tdom,0.50)) ->sumstats_dom
sumstats_dom |> View()

mod1_dom <- lm(tdom~beta*betap*vaxrate/nu,data=expresults_dom,contrasts = list(beta=contr.SAS,
                                                                                  betap=contr.SAS,
                                                                                  nu=contr.SAS))

# mod1_dom <- lm(tdom~beta*betap*vaxrate/nu*mu,data=expresults_dom,contrasts = list(beta=contr.SAS,
#                                                                            betap=contr.SAS,
#                                                                            nu=contr.SAS))


eff_size(emmeans(mod1_dom,pairwise~betap|beta*vaxrate,adjust="tukey"),
         sigma=sigma(mod1_dom),
         edf=df.residual(mod1_dom))


mod1_dom |> aov() |> summary()
emmeans(mod1_dom,~beta*betap*vaxrate/nu)
# emmeans(mod1_dom,~beta*betap*vaxrate/mu)

emmeans(mod1_dom,pairwise~nu|beta*betap*vaxrate,adjust="tukey")
# emmeans(mod1_dom,pairwise~mu|beta*betap*vaxrate,adjust="tukey")
emmeans(mod1_dom,pairwise~betap|beta*vaxrate,adjust="tukey")


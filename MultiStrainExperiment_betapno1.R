library(tidyverse)
library(emmeans)
library(car)
library(gmodels)

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
    mutate(infected = i1+i2+0.000001,infprop1=i1/infected,infprop2=i2/infected) -> cases 
  
  return(cases)
}

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
    tdominance <- NA #never dominates
  }
  else{
    #initialize the logistic model using logit transformed model estimates
    coef(lm(logit(infprop2)~time,data=trun)) |> unname()->tcoeffs
    
    phi_2 <- tcoeffs[1]
    phi_3<- tcoeffs[2]
    
    logisticmodel <- tryCatch({nls(infprop2~SSlogis(time,Asym=1,xmid,scal),
                                   data=trun,
                                   control=list(maxiter=500))},
                              error=function(e){
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


#ORIGINAL STRAIN COULD BE gamma =0.7, gammap = 0.7
gamma <- 0.7
gammap <- 0.7
betaset <- c(1.4,1.75,2.1)
betapset <- c(1.25,1.5,1.75,2) #multiple of beta
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



# rm(expresults)
# sim_array |> rowwise() |> 
#   mutate(tdom = exp_trun(c(SimNumber,beta,betap,gamma,gammap,vaxrate,
#                            ini_recovered,ini_infected,nu,mu))) |>  ungroup()|> 
#   mutate_at(vars(beta,betap,gamma,gammap,vaxrate,
#                  ini_recovered,ini_infected,nu,mu),factor)-> expresults
# 
# 
# write.csv(expresults,file="simresults.csv")

expresults <- read.csv("simresults.csv")

expresults |>  group_by(beta,betap,gamma,gammap,vaxrate,
                        ini_recovered,ini_infected,nu,mu) |>
  mutate(nondom=ifelse(is.na(tdom),1,0)) |> 
  summarise(mean=mean(tdom,na.rm=T),sd=sd(tdom,na.rm=T),max=max(tdom,na.rm=T),min=min(tdom,na.rm=T),
            prop_nondom=sum(nondom,na.rm=T)/n(),median=quantile(tdom,0.50,na.rm=T)) ->sumstats

sumstats  |> View()
mod1 <- lm(tdom~beta*betap*vaxrate*nu,data=expresults,
           contrasts = list(beta=contr.SAS,
                            betap=contr.SAS,
                            nu=contr.SAS))


# mod1 <- lm(tdom~beta*betap*vaxrate*nu*mu,data=expresults,contrasts = list(beta=contr.SAS,
#                                                                betap=contr.SAS,
#                                                                nu=contr.SAS))
mod1 |> aov() |> summary()
emmeans(mod1,~nu*beta*betap*vaxrate)
emmeans(mod1,pairwise~nu|beta*betap*vaxrate,adjust="tukey")


emmeans(mod1,~beta*betap*vaxrate*nu) ->mod1means
mod1means |>  as.data.frame() |> rename(initialvaxx=vaxrate) |>
  ggplot(aes(x=betap,y=emmean, group=nu, shape=nu)) + 
  facet_wrap(~initialvaxx + beta,nrow=2, labeller = label_both) + 
  theme_bw()+
  geom_point(size=2) + geom_line()+ geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL, width=0.1))+
  scale_shape_manual(values=c(15,16))+
  ylab("TTD LSMeans") + xlab("ESTR") + labs(shape="Vaccination Rate") ->lsmeansplot
ggsave(lsmeansplot,filename="lsmeans_no1.png", width=6, height=4, units="in")

# emmeans(mod1_dom,~beta*betap*vaxrate/mu)

# emmeans(mod1_dom,pairwise~nu|beta*betap*vaxrate) |>
#   confint()
# 
# emmeans(mod1_dom,pairwise~betap|beta*nu*vaxrate) |>
#   confint(adjust="bonferroni")






emmeans(mod1,pairwise~nu|betap*vaxrate) |>
  confint(adjust="bonferroni") 

emmeans(mod1,pairwise~betap|nu*vaxrate) |>
  confint(adjust="bonferroni")


emmeans(mod1,poly~betap|nu*vaxrate)

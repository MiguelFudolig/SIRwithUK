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

tdom <- function(ini, x, timestart=0, timeend=10,tol=1e-6){
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
  
  if(tail(trun$infprop2,1) <=tol & max(trun$infprop2) <=0.5){
    tdominance <- NA #never dominates
  }
  else{
    #initialize the logistic model using logit transformed model estimates
    # coef(lm(logit(infprop2)~time,data=trun)) |> unname()->tcoeffs
    # 
    # phi_2 <- tcoeffs[1]
    # phi_3<- tcoeffs[2]
    
    logisticmodel <- tryCatch({nls(infprop2~SSlogis(time,Asym,xmid,scal),
                                   data=filter(trun,infprop2 > tol),
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


#no dominance

initial_recovered <- 0.01 #approximation
initial_infected <- 0.01 
vaxrate <- 0.75
ini <- c(s=1-vaxrate-initial_infected-initial_recovered,
         i1=0.99*initial_infected,
         i2=0.01*initial_infected,
         v=vaxrate,
         r1=initial_recovered,
         r2=0,
         d=0)

x <- list(beta=1.4, 
          betap=1.25*1.4,
          gamma=0.7, 
          gammap=0.7,
          nu=0.25, 
          mu=0)
mod_sir(ini=ini,
        x=x,
        timestart=0,
        timeend = 12) -> results

results |> select(time, infprop1,infprop2) |> 
  rename(Existing=infprop1, Emergent=infprop2) |> 
  pivot_longer(c(Existing,Emergent),names_to="Strain", values_to="Proportion") |> 
  mutate_at(vars(Strain),factor) |> 
  ggplot(aes(x=as.factor(time),y=Proportion,group=Strain, shape=Strain)) +
  theme_bw() + 
  geom_point(size=4) + 
  geom_line(linewidth=1.25) + 
  labs(x="Time", y="Proportion of Infected") ->fail_plot0

fail_plot0 |> ggsave(filename = "failure_0.png",width=6,height=4,units="in")

tdom(ini=ini,x=x,
     timestart=0,
     timeend = 12)


initial_recovered <- 0.01 #approximation
initial_infected <- 0.01 
vaxrate <- 0
ini <- c(s=1-vaxrate-initial_infected-initial_recovered,
         i1=0.99*initial_infected,
         i2=0.01*initial_infected,
         v=vaxrate,
         r1=initial_recovered,
         r2=0,
         d=0)

x <- list(beta=1.4, 
          betap=1.25*1.4,
          gamma=0.7, 
          gammap=0.7,
          nu=0, 
          mu=0)
mod_sir(ini=ini,
        x=x,
        timestart=0,
        timeend = 12) -> results

results |> select(time, infprop1,infprop2) |> 
  rename(Existing=infprop1, Emergent=infprop2) |> 
  pivot_longer(c(Existing,Emergent),names_to="Strain", values_to="Proportion") |> 
  mutate_at(vars(Strain),factor) |> 
  ggplot(aes(x=as.factor(time),y=Proportion,group=Strain, shape=Strain)) +
  theme_bw() + 
  geom_point(size=4) + 
  geom_line(linewidth=1.25) + 
  labs(x="Time", y="Proportion of Infected") ->fail_plot1

fail_plot1 |> ggsave(filename = "failure_1.png",width=6,height=4,units="in")

tdom(ini=ini,x=x,
     timestart=0,
     timeend = 12)

#second way of failure

initial_recovered <- 0.01 #approximation
initial_infected <- 0.01 
vaxrate <- 0.75#0.3
ini <- c(s=1-vaxrate-initial_infected-initial_recovered,
         i1=0.99*initial_infected,
         i2=0.01*initial_infected,
         v=vaxrate,
         r1=initial_recovered,
         r2=0,
         d=0)

x <- list(beta=2.1, 
          betap=2.*2.1,
          gamma=0.7, 
          gammap=0.7,
          nu=0.3, 
          mu=0)
mod_sir(ini=ini,
        x=x,
        timestart=0,
        timeend = 12) -> results2

results2 |> select(time, infprop1,infprop2) |> 
  rename(Existing=infprop1, Emergent=infprop2) |> 
  pivot_longer(c(Existing,Emergent),names_to="Strain", values_to="Proportion") |> 
  mutate_at(vars(Strain),factor) |> 
  ggplot(aes(x=as.factor(time),y=Proportion,group=Strain, shape=Strain)) +
  theme_bw() + 
  geom_point(size=4) + 
  geom_line(linewidth=1.25) + 
  labs(x="Time", y="Proportion of Infected") ->fail_plot2

fail_plot2 |> ggsave(filename = "failure_2.png",width=6,height=4,units="in")

nls(infprop2~SSlogis(time,Asym,xmid,scal),
    data=filter(results2,infprop2 >1e-9),
    control=list(maxiter=500)) ->coeffs

tdom(ini=ini,x=x,
     timestart=0,
     timeend = 12)


#NORMAL LOGISTIC CURVE


initial_recovered <- 0.01 #approximation
initial_infected <- 0.01 
vaxrate <- 0.3
ini <- c(s=1-vaxrate-initial_infected-initial_recovered,
         i1=0.99*initial_infected,
         i2=0.01*initial_infected,
         v=vaxrate,
         r1=initial_recovered,
         r2=0,
         d=0)

x <- list(beta=1.4, 
          betap=1.5*1.4,
          gamma=0.7, 
          gammap=0.7,
          nu=0.3, 
          mu=0)
mod_sir(ini=ini,
        x=x,
        timestart=0,
        timeend = 12) -> results2

results2 |> select(time, infprop1,infprop2) |> 
  rename(Existing=infprop1, Emergent=infprop2) |> 
  pivot_longer(c(Existing,Emergent),names_to="Strain", values_to="Proportion") |> 
  mutate_at(vars(Strain),factor) |> 
  ggplot(aes(x=as.factor(time),y=Proportion,group=Strain, shape=Strain)) +
  theme_bw() + 
  geom_point(size=4) + 
  geom_line(linewidth=1.25) + 
  labs(x="Time", y="Proportion of Infected") ->okay_plot2

okay_plot2 |> ggsave(filename = "okay.png",width=6,height=4,units="in")


tdom(ini=ini,x=x,
     timestart=0,
     timeend = 12)


#ORIGINAL STRAIN COULD BE gamma =0.7, gammap = 0.7
gamma <- 0.7
gammap <- 0.7
betaset <- c(1.4,1.75,2.1)
betapset <- c(1.25,1.5,1.75,2) #multiple of beta
initial_recovered <- 0.01 #approximation
initial_infected <- 0.01 #loosely based on data
vaxxrateset <- c(0,0.25, 0.5, 0.75) #c(0,0.3)
nuset <- c(0,0.25, 0.5) #c(0, 0.3)
muset <- c(0)
nsims <-10 #100

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


#uncomment to simulate results
rm(expresults)
sim_array |> rowwise() |>
  mutate(tdom = exp_trun(c(SimNumber,beta,betap,gamma,gammap,vaxrate,
                           ini_recovered,ini_infected,nu,mu))) |>  ungroup()|>
  mutate_at(vars(beta,betap,gamma,gammap,vaxrate,
                 ini_recovered,ini_infected,nu,mu),factor)-> expresults
expresults |> filter(is.na(tdom))

write.csv(expresults,file="simresults_afterrevision_3X4X4X3.csv")

expresults <- read.csv("simresults_afterrevision_3X4X4X3.csv")

expresults |> 
    mutate_at(vars(beta,betap,gamma,gammap,vaxrate,
                   ini_recovered,ini_infected,nu,mu),factor)-> expresults
  
expresults |>  group_by(beta,betap,gamma,gammap,vaxrate,
                        ini_recovered,ini_infected,nu,mu) |>
  mutate(nondom=ifelse(is.na(tdom),1,0)) |> 
  summarise(mean=mean(tdom,na.rm=T),sd=sd(tdom,na.rm=T),max=max(tdom,na.rm=T),min=min(tdom,na.rm=T),
            prop_nondom=sum(nondom,na.rm=T)/n(),median=quantile(tdom,0.50,na.rm=T)) ->sumstats

sumstats  |> View()

sumstats |> ungroup() |> filter(prop_nondom > 0) |>
  select(beta,betap,nu,vaxrate,prop_nondom) |> View()






mod1 <- lm(tdom~beta*betap*vaxrate*nu,data=expresults,
           contrasts = list(beta=contr.SAS,
                            betap=contr.SAS,
                            nu=contr.SAS,
                            vaxrate=contr.SAS))

#only second order
mod2 <- lm(tdom~beta+betap+vaxrate+nu +
             beta:betap + beta:vaxrate + beta:nu + 
             betap:vaxrate + betap:nu + vaxrate:nu,data=expresults,
                   contrasts = list(beta=contr.SAS,
                                    betap=contr.SAS,
                                    nu=contr.SAS,
                                    vaxrate=contr.SAS))


# mod1 <- lm(tdom~beta*betap*vaxrate*nu*mu,data=expresults,contrasts = list(beta=contr.SAS,
#                                                                betap=contr.SAS,
#                                                                nu=contr.SAS))
plot(mod1)
#full model
mod1 |> aov() |> summary()
#second order interaction only
mod2 |> aov() |> summary()


emmeans(mod1,~nu*beta*betap*vaxrate)
emmeans(mod1,pairwise~nu|beta*betap*vaxrate,adjust="tukey")


emmeans(mod1,~beta*betap*vaxrate*nu) ->mod1means
mod1means |>  as.data.frame() |> rename(initialvaxx=vaxrate) |>
  ggplot(aes(x=betap,y=emmean, group=nu, shape=nu)) + 
  facet_grid(beta~initialvaxx, labeller = label_both) + 
  theme_bw()+
  geom_point(size=3) + geom_line()+ geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL, width=0.3))+
  scale_shape_manual(values=c(15,16,17))+
  ylab("TTD LSMeans") + xlab("ESTR") + labs(shape="Vaccination Rate") ->lsmeansplot
ggsave(lsmeansplot,filename="./3X4X4X3/lsmeans_no1.png", width=12, height=9, units="in")

# emmeans(mod1_dom,~beta*betap*vaxrate/mu)

# emmeans(mod1_dom,pairwise~nu|beta*betap*vaxrate) |>
#   confint()
# 
# emmeans(mod1_dom,pairwise~betap|beta*nu*vaxrate) |>
#   confint(adjust="bonferroni")




#vaccination rate

emmeans(mod1,pairwise~nu|betap*vaxrate*beta) |>
  confint(adjust="bonferroni") ->nusimp

nusimp |> as.data.frame() |> 
  write.csv("./3X4X4X3/nu_simpleeffects.csv")


emmeans(mod1,pairwise~nu|betap*vaxrate*beta)$contrasts |> 
  confint(adjust="bonferroni") |>
  as.data.frame() |> rename(initialvaxx=vaxrate) |> 
  ggplot(aes(x=betap,y=estimate, group=beta, shape=beta)) + 
  theme_bw() + 
  facet_wrap(~initialvaxx, labeller = label_both) + 
  geom_point(size=2) + geom_line()+ geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL, width=0.3))+
  scale_shape_manual(values=c(15,16,17))+
  ylab("TTD Difference") + xlab("ESTR") + labs(shape="Existing Transmission") ->vaxratediff

vaxratediff |> ggsave(filename="./3X4X4X3/vaxrate_diffplot.png",width=6, height=4, units="in")



#initial vaccination

emmeans(mod1,pairwise~vaxrate|nu*betap*beta) |>
  confint(adjust="bonferroni")-> vaxsimp

vaxsimp$contrasts|> as.data.frame() |> 
  write.csv("./3X4X4X3/vaxrate_simpleeffects.csv")

emmeans(mod1,pairwise~vaxrate|nu*betap*beta)$contrasts |> 
  confint(adjust="bonferroni") |>
  as.data.frame() |> 
  ggplot(aes(x=betap,y=estimate, group=beta, shape=beta)) + 
  theme_bw() + 
  facet_wrap(~nu, labeller = label_both) + 
  geom_point(size=2) + geom_line()+ geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL, width=0.3))+
  scale_shape_manual(values=c(15,16,17))+
  ylab("TTD Difference") + xlab("ESTR") + labs(shape="Existing Transmission") ->initialvaxdiff


initialvaxdiff |> ggsave(filename="./3X4X4X3/initialvaxrate_diffplot.png",width=6, height=4, units="in")


emmeans(mod1,pairwise~vaxrate|nu*betap*beta) |>
  confint(adjust="bonferroni") |> as.data.frame() |> 
  summarise(min=min(contrasts.estimate),max=max(contrasts.estimate))

#ESTR

joint_tests(mod1,by=c("nu","vaxrate","beta")) |> write.csv("./3X4X4X3/estr_sliceeffects.csv")

#two-way interaction between beta and beta'
joint_tests(mod1,by=c("nu","vaxrate")) |> write.csv("./3X4X4X3/estr_interaction.csv")

joint_tests(mod1,by=c("beta")) |> write.csv("./3X4X4X3/estr_three-wayinteraction.csv")


emmeans(mod1,pairwise~betap|nu*vaxrate*beta) |>
  confint(adjust="bonferroni") -> estrcontrast

estrcontrast$contrasts |> as.data.frame() |> 
  write.csv("./3X4X4X3/estr_simpleeffects.csv")

estrcontrast$emmeans |> as.data.frame() |> 
  write.csv("./3X4X4X3/estr_lsmeans.csv")


estrcontrast$contrasts |>
  as.data.frame() |> rename(initialvaxx=vaxrate) |> 
  ggplot(aes(x=beta,y=estimate, group=nu, shape=nu,color=nu)) + 
  theme_bw() + 
  facet_wrap(contrast~initialvaxx, labeller = label_both) + 
  geom_point(size=3) + geom_line()+ geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL, width=0.5))+
  scale_shape_manual(values=c(15,16))+
  ylab("TTD Difference") + xlab("Existing Strain Transmission Coefficient") + 
  labs(shape="Vaccination Rate", color="Vaccination Rate") ->estrdiffplot
estrdiffplot |> 
  ggsave(filename="./3X4X4X3/estr_diffplot.png",width=12,height=12,units="in")

estrcontrast$contrasts |>
   as.data.frame() |> 
  filter(vaxrate==0.75 | nu==0.5) |> summarise(min=min(estimate),max=max(estimate))

emmeans(mod1,poly~betap|nu*vaxrate*beta) |> as.data.frame() |> 
  write.csv("./3X4X4X3/estr_quadratic.csv")


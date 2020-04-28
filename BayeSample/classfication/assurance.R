

######################################################################### 
# 1.1. Use R "integrate" instead of the trapezoidal numerical integration 
#   in  Chuang-Stein (2006) which can compute the assurance faster 
######################################################################### 




sprob.Chuang = function(prior.mean,prior.sd,prior.size,post.size){ 
  prior.sdm = sqrt(2/prior.size)*prior.sd # prior sd for the mean 
  post.sdm  = sqrt(2/post.size)*prior.sd # posterior sd for the mean 
  # fn for the Prob(trial produces a significant p-val)*prior distribution 
  integrand <- function(delta) 
    pnorm(1.96*post.sdm,mean=delta,sd=post.sdm,lower.tail=FALSE,log.p=FALSE)* 
    dnorm(delta,mean=prior.mean,sd=prior.sdm,log = FALSE) 
  # Numerical integration of delta from -Inf to Inf 
  avg = integrate(integrand, lower = -Inf, upper = Inf)$value 
  # output 
  avg 
}

inv.sprob.Chuang = function(prior.mean,prior.sd,prior.size,assurance,Textparam=NULL,CNname=c("先验均值","先验方差","先验样本","assurance")){ 
  
  #prior.mean <- 2.5
  #prior.sd <-7.14
  #prior.size<-25
  #assurance <- 0.8
  u <- function(post.size){
    prior.sdm = sqrt(2/prior.size)*prior.sd # prior sd for the mean 
    post.sdm  = sqrt(2/post.size)*prior.sd # posterior sd for the mean 
    # fn for the Prob(trial produces a significant p-val)*prior distribution 
    integrand <- function(delta) 
      pnorm(1.96*post.sdm,mean=delta,sd=post.sdm,lower.tail=FALSE,log.p=FALSE)* 
      dnorm(delta,mean=prior.mean,sd=prior.sdm,log = FALSE) 
    # Numerical integration of delta from -Inf to Inf 
    avg = integrate(integrand, lower = -Inf, upper = Inf)$value 
  }
  #最终结果的范围必须再lower-upper之间
  #asurace总是小于
  sqrt.inv = inversed(u, lower=0, upper=30000)
  sample <- sqrt.inv(assurance)   
  new.sample <- list(prior.mean=prior.mean,prior.sd=prior.sd,prior.size=prior.size,assurance=assurance, sample=sample)
}


library(GoFKernel)
u <- function(post.size){
  prior.mean <-2.5
  prior.sd <-7.14
  prior.size <- 25
  prior.sdm = sqrt(2/prior.size)*prior.sd # prior sd for the mean 
  post.sdm  = sqrt(2/post.size)*prior.sd # posterior sd for the mean 
  # fn for the Prob(trial produces a significant p-val)*prior distribution 
  integrand <- function(delta) 
    pnorm(1.96*post.sdm,mean=delta,sd=post.sdm,lower.tail=FALSE,log.p=FALSE)* 
    dnorm(delta,mean=prior.mean,sd=prior.sdm,log = FALSE) 
  # Numerical integration of delta from -Inf to Inf 
  avg = integrate(integrand, lower = -Inf, upper = Inf)$value 
}

#a  = u(127)
#f.inv <- inverse(u,lower=-Inf,upper=Inf)
#f.inv(0.63)


inversed = function(fn, interval = NULL, lower = min(interval), upper = max(interval), ...){
  Vectorize(function(y){
    uniroot(f=function(x){fn(x)-y}, lower=lower, upper=upper, ...)$root
  })
}

u <- function(post.size){
  prior.mean <-2.5
  prior.sd <-0.1
  prior.size <- 25
  prior.sdm = sqrt(2/prior.size)*prior.sd # prior sd for the mean 
  post.sdm  = sqrt(2/post.size)*prior.sd # posterior sd for the mean 
  # fn for the Prob(trial produces a significant p-val)*prior distribution 
  integrand <- function(delta) 
    pnorm(1.96*post.sdm,mean=delta,sd=post.sdm,lower.tail=FALSE,log.p=FALSE)* 
    dnorm(delta,mean=prior.mean,sd=prior.sdm,log = FALSE) 
  # Numerical integration of delta from -Inf to Inf 
  avg = integrate(integrand, lower = -Inf, upper = Inf)$value 
}

#x = 127
#f <- function(x){
#  x**2
#}
#y = 0.8
#最终结果的范围必须再lower-upper之间
#asurace总是小于
#sqrt.inv = inversed(u, lower=0, upper=30000)
#sqrt.inv(y)



##################################################################### 
# 2.1: Assurance for Normal Data with Known Sigma (ANDks)  
#  to check with the results from Appendix 1 
###################################################################### 
ANDks = function(nsimu=10000,prior.mean,prior.sd,prior.size,post.size,alpha){ 
  sim.pow = rep(0, nsimu) 
  for(i in 1:nsimu){ 
    # calculate the standard deviation for the means 
    prior.sdm = sqrt(2/prior.size)*prior.sd # prior sd for the mean 
    post.sdm  = sqrt(2/post.size)*prior.sd # posterior sd for the mean 
    # sample the prior 
    Delta = rnorm(1,prior.mean,prior.sdm) 
    # with the sampled prior, calculate the power 
    sim.pow[i]=pnorm(qnorm(1-alpha/2)*post.sdm, 
                     mean=Delta,sd=post.sdm,lower.tail=FALSE,log.p=FALSE) 
  } # end of i-loop 
  # average the simulated power 
  mean(sim.pow) 
} # end of "ANDks" function 
## run the code to check with the calculations in Appendix 1 




inv.ANDks <- function(prior_mean,prior_sd,prior_size,N,alpha=0.05, Textparam=NULL,CNname=c("先验均值","先验方差","先验样本","样本量","α"),
                      tips=list(prior_mean="正态分布的第一个参数,假设两组疗效差异服从该正态分布",
                                prior_sd="正态分布的第二个参数,假设两组疗效差异服从该正态分布",
                                prior_size="先验信息能提供的样本量",
                                N="传统方法所需的样本量",
                                alpha="I类错误的概率"
                      ),
                      defaultValue=list(prior_mean=2.5,prior_sd=7.14,prior_size=25,alpha=0.05)
){
  
  nsimu <- 10000
  u <- function(post.size){
    sim.pow = rep(0, nsimu) 
    for(i in 1:nsimu){ 
      # calculate the standard deviation for the means 
      prior_sdm = sqrt(2/prior_size)*prior_sd # prior sd for the mean 
      post.sdm  = sqrt(2/post.size)*prior_sd # posterior sd for the mean 
      # sample the prior 
      Delta = rnorm(1,prior_mean,prior_sdm) 
      # with the sampled prior, calculate the power 
      sim.pow[i]=pnorm(qnorm(1-alpha/2)*post.sdm, 
                       mean=Delta,sd=post.sdm,lower.tail=FALSE,log.p=FALSE) 
    } # end of i-loop 
    # average the simulated power 
    mean(sim.pow)   
  }
  
  #由assurance反函数求N
  #此时给定的参数为assurance而不是post.size
  # sqrt.inv = inversed(u, lower=0, upper=30000)
  # sample <- sqrt.inv(assurance)   
  assurance <- round(mapply(u, N),3)
  
  conclusion <- paste0("在α为",alpha,",样本量为",N[1],"，两组均数差异服从N(",prior_mean,"，",prior_sd,")，先验样本量为",prior.size,"时，试验的平均功效为",assurance[1])
  new.sample <- list(prior_mean=prior_mean,prior_sd=prior_sd,assurance=assurance, alpha=alpha,sample=N,conclusion=conclusion)
  
}

#inv.ANDks(1000,2.5,7.14,25,c(0.6,0.7,0.8))


#invANDks(1000,2.5,7.14,25,c(0.6,0.7,0.8,0.85))
###################################################################### 
# 2.2: Assurance for Normal Data with unknown sigma(ANDus) 
###################################################################### 
ANDus = function(nsimu,prior.mean,prior.sd,prior.size,post.size){ 
  sim.pow = rep(0, nsimu) 
  for(i in 1:nsimu){ 
    # sample chisq for sigma since (n-1)*s^2/sigma^2 ~chisq(n-1) 
    sd = sqrt((prior.size-1)*prior.sd^2/rchisq(1,df=prior.size-1)) 
    # calculate the standard deviation for the mean 
    prior.sdm = sqrt(2/prior.size)*sd # prior sd for the mean 
    post.sdm  = sqrt(2/post.size)*sd # posterior sd for the mean 
    # sample the prior 
    Delta = rnorm(1, prior.mean,prior.sdm) 
    # with the sampled prior, calculate the power 
    sim.pow[i] = pnorm(1.96*post.sdm,mean=Delta, sd=post.sdm,lower.tail=FALSE,log.p=FALSE) 
  } # end of i-loop 
  # assurance is the average of simulated power 
  mean(sim.pow) 
} # end of "ANDus" function 


inv.ANDus = function(prior_mean,prior_sd,prior_size,N,alpha,Textparam=NULL,CNname=c("先验均值","先验方差","先验样本","样本量","α"),
                     tips=list(prior_mean="正态分布的第一个参数,假设两组疗效差异服从该正态分布",
                               prior_sd="正态分布的第二个参数,假设两组疗效差异服从该正态分布",
                               prior_size="先验信息能提供的样本量",
                               N="传统方法所需的样本量",
                               alpha="I类错误的概率"
                     ),
                     defaultValue=list(prior_mean=2.5,prior_sd=7.14,prior_size=25,alpha=0.05)
                     ){ 
  nsimu <- 10000
  u <- function(post.size){
    sim.pow = rep(0, nsimu) 
    for(i in 1:nsimu){ 
      # sample chisq for sigma since (n-1)*s^2/sigma^2 ~chisq(n-1) 
      sd = sqrt((prior_size-1)*prior_sd^2/rchisq(1,df=prior_size-1)) 
      # calculate the standard deviation for the mean 
      prior_sdm = sqrt(2/prior_size)*sd # prior sd for the mean 
      post.sdm  = sqrt(2/post.size)*sd # posterior sd for the mean 
      # sample the prior 
      Delta = rnorm(1, prior_mean,prior_sdm) 
      # with the sampled prior, calculate the power 
      sim.pow[i] = pnorm(qnorm(1-alpha/2)*post.sdm,mean=Delta, sd=post.sdm,lower.tail=FALSE,log.p=FALSE) 
    } # end of i-loop 
    # assurance is the average of simulated power 
    mean(sim.pow) 
  }
  
  #sqrt.inv = inversed(u, lower=0, upper=30000)
  #sample <- sqrt.inv(assurance)   
  assurance <- round(mapply(u, N),3)
  
  conclusion <- paste0("在α为",alpha,",样本量为",N[1],"，两组均数差异服从N(",prior_mean,"，",prior_sd,")，先验样本量为",prior.size,"时，试验的平均功效为",assurance[1])
  new.sample <- list(prior_mean=prior_mean,prior_sd=prior_sd,prior_size=prior_size,assurance=assurance, sample=N,conclusion=conclusion)
  
} 



#a = inv.ANDus(1000,2.5,7.14,25,c(0.6,0.7,0.8))
######################################################################### 
# 2.3: Assurance with binary clinical trial in O‟Hagan et al. Example 4 
########################################################################## 
######################################################################### 
# 2.3: Assurance with binary clinical trial in O’Hagan et al. Example 4 
########################################################################## 
AWBD <- function(nsimu=10000,alpha=0.05,k=1, beta1param1, beta1param2,beta2param1, beta2param2,post.size1,type="RD",Textparam=NULL,CNname=c("模拟次数","α","分配比例","a/先验1-beta(a,b)","b/先验1-beta(a,b)","a/先验2-beta(a,b)","b/先验2-beta(a,b)","后验样本","比较类型")){
    sim.pow = rep(0,nsimu) 
    for (i in 1:nsimu){ 
      p1 = rbeta(1,beta1param1,beta1param2) 
      # sample Beta for p1  p1 = rbeta(1,5,20) 
      # sample Beta for p2 from a mixture  
      p2 = rbeta(1,beta2param1,beta2param2) 
      # z-value in equation  
      if(type=="RD"){
        z.val = (p2-p1)/sqrt(p1*(1-p1)/post.size1+p2*(1-p2)/(post.size1/k)) 
      }else{
        OR = (p2/(1-p2))/(p1/(1-p1))
        z.val = log(OR)/sqrt(1/(post.size1*p1*(1-p1))+1/(post.size1*p2*(1-p2)))
      }
      # with the sampled prior, calculate the power 
      sim.pow[i] = pnorm(-qnorm(1-alpha/2)+z.val) 
    } # end of i-loop 
    # Assurance as the mean 
    mean(sim.pow) 
  }


inv.AWBD <- function(beta1param1, beta1param2,beta2param1, beta2param2,N,alpha=0.05,Textparam=NULL,CNname=c("a/先验1-beta(a,b)","b/先验1-beta(a,b)","a/先验2-beta(a,b)","b/先验2-beta(a,b)","样本量","α"),
                     tips=list(beta1param1="beta分布的第一个参数，组1率服从给beta分布",
                               beta1param2="beta分布的第一个参数，组1率服从给beta分布",
                               beta2param1="beta分布的第一个参数，组2率服从给beta分布",
                               beta2param2="beta分布的第一个参数，组2率服从给beta分布",
                               N="传统方法所需的样本量",
                               alpha="I类错误的概率"
                     ),
                     defaultValue=list(beta1param1=1, beta1param2=1,beta2param1=1, beta2param2=2,alpha=0.05)
                     ){
  nsimu=10000
  k=1
  type="OR"
  u <- function(post.size1){
    sim.pow = rep(0,nsimu) 
    for (i in 1:nsimu){ 
      p1 = rbeta(1,beta1param1,beta1param2) 
      # sample Beta for p1  p1 = rbeta(1,5,20) 
      # sample Beta for p2 from a mixture  
      p2 = rbeta(1,beta2param1,beta2param2) 
      # z-value in equation  
      if(type=="RD"){
        z.val = (p2-p1)/sqrt(p1*(1-p1)/post.size1+p2*(1-p2)/(post.size1/k)) 
      }else{
        OR = (p2/(1-p2))/(p1/(1-p1))
        z.val = log(OR)/sqrt(1/(post.size1*p1*(1-p1))+1/(post.size1*p2*(1-p2)))
      }
      # with the sampled prior, calculate the power 
      sim.pow[i] = pnorm(-qnorm(1-alpha/2)+z.val) 
    }
    mean(sim.pow) 
  }
  
  #sqrt.inv = inversed(u, lower=0, upper=30000)
  #sample <- sqrt.inv(assurance)   
  assurance <- round(mapply(u, N),3)
  conclusion <- paste0("在α为",alpha,",样本量为",N[1],"，两组率分布服从p1~beta(",beta1param1,"，",beta1param2,"),p2~beta(",beta2param1,",",beta2param2,")时，试验的平均功效为",assurance[1])
  new.sample <- list(alpha=alpha,beta1param1=beta1param1, beta1param2=beta1param2,beta2param1=beta2param1, beta2param2=beta2param2,assurance=assurance, sample=N,conclusion=conclusion)
  
}

#a = inv.AWBD(nsimu=1000,alpha=0.05,k=0.7,beta1param1=5, beta1param2=10,beta2param1=3, beta2param2=4.5, p1=0.23,p2=0.12,assurance=c(0.5,0.55))

library(Rlab) # for rbern 
nsimu = 1000000;alpha=0.05; post.size1 = 200; post.size2 = 400 
sim.pow = rep(0,nsimu) 
for (i in 1:nsimu){ 
  # sample Beta for p1  p1 = rbeta(1,5,20) 
  # sample Beta for p2 from a mixture  
  w = rbern(1,0.15);p2=w*rbeta(1,2,23)+(1-w)*rbeta(1,3,4.5) 
  # z-value in equation  
  z.val = (p2-p1)/sqrt(p1*(1-p1)/post.size1+p2*(1-p2)/post.size2) 
  # with the sampled prior, calculate the power 
  sim.pow[i] = pnorm(-qnorm(1-alpha/2)+z.val) 
} # end of i-loop 
# Assurance as the mean 
mean(sim.pow) 

######################################################## 
# MCSB function  
######################################################## 
pow2assurance=function(nsimu,n1,n2,pA,pB,sig2A,sig2B,sig2A2,sig2B2,alpha){ 
  # initializes the power and assurance 
  pow = assu1 = assu2 = 0 
  # loop for calculation 
  for(i in 1:nsimu){ 
    ### power simulation 
    xA=rbinom(n1,1,pA);xB=rbinom(n2,1,pB) 
    dd=data.frame(x=c(xA,xB),trt=c(rep("A",n1),rep("B",n2))) 
    md=glm(x~trt,dd,family="binomial"); 
    pval.md = summary(md)$coef["trtB","Pr(>|z|)"] 
    pow = pow+sum(pval.md < alpha) 
    ### assurance simulation for sigma1 
    alphaA=rnorm(1,log(pA/(1-pA)),sqrt(sig2A)); 
    alphaB=rnorm(1,log(pB/(1-pB)),sqrt(sig2B)) 
    pAs = exp(alphaA)/(1+exp(alphaA));pBs = exp(alphaB)/(1+exp(alphaB)); 
    xA = rbinom(n1,1,pAs);xB = rbinom(n2,1,pBs) 
    dd = data.frame(x = c(xA, xB), trt = c(rep("A", n1), rep("B", n2))) 
    md = glm(x~trt, dd, family="binomial"); 
    pval.md = summary(md)$coef["trtB","Pr(>|z|)"] 
    assu1 = assu1+sum(pval.md < alpha) 
    ### assurance simulation for sigma2 
    alphaA=rnorm(1,log(pA/(1-pA)),sqrt(sig2A2)); 
    alphaB=rnorm(1,log(pB/(1-pB)),sqrt(sig2B2)) 
    pAs = exp(alphaA)/(1+exp(alphaA));pBs = exp(alphaB)/(1+exp(alphaB)); 
    xA = rbinom(n1,1,pAs);xB = rbinom(n2,1,pBs) 
    dd = data.frame(x = c(xA, xB), trt = c(rep("A", n1), rep("B", n2))) 
    md = glm(x~trt, dd, family="binomial"); 
    pval.md = summary(md)$coef["trtB","Pr(>|z|)"] 
    assu2 = assu2+sum(pval.md < alpha) 
  } #End of i-loop 
  # output 
  list(pow = pow/nsimu,assu1=assu1/nsimu,assu2=assu2/nsimu ) 
} # end of pow2assurance 


######################################################################### 
# 4. MCSB:Exponential distribution  
######################################################################### 
#a = assurance.exp(betaparameters=c(1,2),normalparameters=c(1,3),t0=2,Q1=0.5,R.trial=3,T.trial=5,n=100)
###### Probability of the event in the power and assurnace calculation #####
power.exponential<-function(s1,t0,theta,R,T,N,Q1){
  lam1<--log(s1)/t0
  lam2<-exp(theta)*lam1
  Pevent<-function(lamda,R,T){
    1-(exp(-lamda*(T-R))-exp(-lamda*T))/(lamda*R)
  }
  
  power<-pnorm(-log(lam2/lam1)/sqrt(1/(N*Q1*Pevent(lam1,R,T))+1/(N*(1-Q1)*Pevent(lam2,R,T)))-qnorm(1-0.025))
  power
}


#######################
### assurance formula #
#######################

### Priors for assurance ###
### s1~Beta(a,b)
### rho~Normal(m,v)

assurance.exponential<-function(t0,a,b,m,v,M,R,T,N,Q1){
  
  rho<-rnorm(M,m,sqrt(v))
  
  s1<-rbeta(M,a,b)
  lam1<--log(s1)/t0
  
  s2<-s1+rho
  pos1<-which(s2<1)
  s2<-s2[pos1]
  pos2<-which(s2>0)
  lam2<--log(s2[pos2])/t0
  lam1<-lam1[pos1][pos2]
  rho<-rho[pos1][pos2]
  
  Pevent<-function(lamda,R,T){
    1-(exp(-lamda*(T-R))-exp(-lamda*T))/(lamda*R)
  }
  
  power<-pnorm(-log(lam2/lam1)/sqrt(1/(N*Q1*Pevent(lam1,R,T))+1/(N*(1-Q1)*Pevent(lam2,R,T)))-qnorm(1-0.025))
  assurance<-sum(power)/M
  assurance
}


inv.ASDExp <- function(t0,a,b,m,v,R,TT,N,alpha,Textparam=NULL,CNname=c("时刻t0","a/先验beta","b/先验beta","先验ρ均值","先验ρ方差","试验组总时间","对照组总时间","样本量","α"),
                       tips=list(t0="用于计算特定时刻的生存率与生存风险的时刻",
                                 a="beta分布的第一个参数，组1率服从给beta分布",
                                 b="beta分布的第一个参数，组1率服从给beta分布",
                                 m="正态分布的第一个参数，对数风险比ρ(ρ=log(λ2/λ1))的先验服从该分布",
                                 v="正态分布的第二个参数，对数风险比ρ(ρ=log(λ2/λ1)的先验服从该分布",
                                 R="试验的招募时间",
                                 TT="试验的总时间长度",
                                 N="样本量",
                                 alpha="I类错误的概率"
                       ),
                       defaultValue=list(t0=5,a=60,b=40,m=0.2,v=0.0001,R=3,TT=5,alpha=0.05)
                       ){
  
  M=10000
  Q1=0.5
  u <- function(N){
    rho<-rnorm(M,m,sqrt(v))
    
    s1<-rbeta(M,a,b)
    lam1<--log(s1)/t0
    
    s2<-s1+rho
    pos1<-which(s2<1)
    s2<-s2[pos1]
    pos2<-which(s2>0)
    lam2<--log(s2[pos2])/t0
    lam1<-lam1[pos1][pos2]
    rho<-rho[pos1][pos2]
    
    Pevent<-function(lamda,R,TT){
      1-(exp(-lamda*(TT-R))-exp(-lamda*TT))/(lamda*R)
    }
    
    power<-pnorm(-log(lam2/lam1)/sqrt(1/(N*Q1*Pevent(lam1,R,TT))+1/(N*(1-Q1)*Pevent(lam2,R,TT)))-qnorm(1-alpha/2))
    assurance<-sum(power)/M
    assurance
  }
  
  #a = u(30000)
  #sqrt.inv = inversed(u, lower=0, upper=30000)
  #sample <- sqrt.inv(assurance)   
  assurance <- round(mapply(u, N),3)
  conclusion <- paste0("在α为",alpha,",样本量为",N[1],"，试验的招募期为",R,"年，随访期为",TT,"年，时刻t0为",t0,"，组1的存活率S1~beta(",a,",",b,"),对数风险比ρ~N(",m,",",v,")时，试验的平均功效为",assurance[1])
  
  new.sample <- list(t0=t0,a=a,b=b,m=m,v=v,R=R,TT=TT,Q1=Q1,assurance=assurance, sample=N,conclusion=conclusion)
  
}

#a = inv.ASDExp(12,1,3,2.5,7.12,1000,3,5,0.05,c(0.05,0.07))


##################################################################################
### Weibull model power and assurance calculation                              ###
### Note:                                                                      ###
### 1. Two independent groups with different                                   ###
###    shape parameters                                                        ###
### 2. Use Gross & Clark (1975)'s power formula                                ###
### 3. Assurance is calculated given priors                                    ###
###   for s11,d11,d12,d21                                                      ###
##################################################################################

### N: the total number of patients in the trial
### Q1: the proportion of patients allocated in the control group
### M: the number of Monte Carlo simulations

### Parameters used to calculate power:
### s11: the survival rate at time t1 in the control group
### s12: the survival rate at time t2 in the control group
### s21: the survival rate at time t1 in the experimental group
### s22: the survival rate at time t2 in the experimental group

### Priors used in the assurance calculation:
### Assume the survival rate at time t1 in the control group has a beta
### distribution, S11~beta(a.s11,b.s11)
### Assume the survival different between time t1 and t2 in the control group
### has a beta distribution d11~beta(a.d11,b.d11)
### Assume the treatment difference between two groups at time t1 has a normal
### distribution, d12~N(m.d12,v.d12)
### Assume the survival different between time t1 and t2 in the experimental group
### has a beta distribution d22~beta(a.d22,b.d22)

assurance.weibull<-function(t1,t2,N,s11,s12,s21,s22,Q1,M,a.s11,b.s11,a.d11,b.d11,
                            m.d12,v.d12,a.d22,b.d22){
  
  ##### Functions for estimating the model parameters #####
  e.gamma<-function(s1,s2,t1,t2){
    log(log(s2)/log(s1))/log(t2/t1)
  }
  
  e.lambda<-function(s1,s2,t1,t2){
    -log(s1)/t1^(e.gamma(s1,s2,t1,t2))
  }
  
  e.mu<-function(s1,s2,t1,t2){
    e.lambda(s1,s2,t1,t2)^(-1/e.gamma(s1,s2,t1,t2))*gamma(1+1/e.gamma(s1,s2,t1,t2))
  }
  
  e.sigma2<-function(s1,s2,t1,t2){
    e.lambda(s1,s2,t1,t2)^(-2/e.gamma(s1,s2,t1,t2))*(gamma(1+2/e.gamma(s1,s2,t1,t2))
                                                     -(gamma(1+1/e.gamma(s1,s2,t1,t2)))^2)
  }
  
  
  ### Power function
  powerW2<-function(s11,s12,s21,s22,t1,t2,N,Q1){
    e.mu1<-e.mu(s11,s12,t1,t2)
    e.mu2<-e.mu(s21,s22,t1,t2)
    e.sigma21<-e.sigma2(s11,s12,t1,t2)
    e.sigma22<-e.sigma2(s21,s22,t1,t2)
    v<-e.sigma21/(N*Q1)+e.sigma22/(N*(1-Q1))
    
    powerW2<-pnorm((e.mu2-e.mu1)/sqrt(v)-qnorm(1-0.05/2))
    powerW2
  }
  
  ### Sample s11,s12,s21,s22 for assurance
  
  sampled.s<-function(M,a.s11,b.s11,a.d11,b.d11,m.d12,v.d12,a.d22,b.d22,t1,t2){
    Smatrix<-matrix(NA,nrow=4,ncol=M)
    for(i in 1:M){
      s12<--1;s21<--1;s22<--1;g1<-10;g2<-10
      while(s12<0||s12>1||s21<0||s21>1||s22<0||s22>1||g1>170||g2>170){
        d11<-rbeta(1,a.d11,b.d11)
        d12<-rnorm(1,m.d12,sqrt(v.d12))
        d22<-rbeta(1,a.d22,b.d22)
        s11<-rbeta(1,a.s11,b.s11)
        
        s12<-s11-d11
        s21<-s11+d12
        s22<-s21-d22
        if(s11<0||s12<0){
          g1<-10
        }else{
          g1<-2/e.gamma(s11,s12,t1,t2)
        }
        
        if(s21<0||s22<0){
          g2<-10
        }else{
          g2<-2/e.gamma(s21,s22,t1,t2)
        }
      }
      
      Smatrix[1,i]<-s11
      Smatrix[2,i]<-s12
      Smatrix[3,i]<-s21
      Smatrix[4,i]<-s22
    }
    Smatrix
  }
  
  
  
  ### Calculate power and assuance
  p<-powerW2(s11,s12,s21,s22,t1,t2,N,Q1)
  sampled.surv<-sampled.s(M,a.s11,b.s11,a.d11,b.d11,m.d12,v.d12,a.d22,b.d22,t1,t2)
  
  ### Simulate the Weibull parameters for two groups
  gamma1<-e.gamma(sampled.surv[1,],sampled.surv[2,],t1,t2)
  lambda1<-e.lambda(sampled.surv[1,],sampled.surv[2,],t1,t2)
  gamma2<-e.gamma(sampled.surv[3,],sampled.surv[4,],t1,t2)
  lambda2<-e.lambda(sampled.surv[3,],sampled.surv[4,],t1,t2)
  
  ### Simulate trial data for each group for each simulated pair of Weibull parameters
  assurance<-function(n,x1,x2,v1,v2){
    x1<-sapply(seq(1:M),function(i){rweibull(n*Q1,shape=gamma1[i],scale=1/lambda1[i])})
    x2<-sapply(seq(1:M),function(i){rweibull(n*(1-Q1),shape=gamma2[i],scale=1/lambda2[i])})
    
    v1<-sapply(seq(1:M),function(i){var(x1[,i])})
    v2<-sapply(seq(1:M),function(i){var(x2[,i])})
    
    test.stats<-(colSums(x2)/dim(x1)[1]-colSums(x1)/dim(x1)[1])/sqrt(v1/n/Q1+v2/n/(1-Q1))
    e.df<-(v1/(n*Q1)+v2/(n*(1-Q1)))^2/(v1^2/(n*Q1)^2/(n*Q1-1)+v2^2/(n*(1-Q1))^2/(n*(1-Q1)-1))
    
    
    while(sum(e.df==0)>0||sum(is.na(e.df))>0){
      x1<-sapply(seq(1:M),function(i){rweibull(n*Q1,shape=gamma1[i],scale=1/lambda1[i])})
      x2<-sapply(seq(1:M),function(i){rweibull(n*(1-Q1),shape=gamma2[i],scale=1/lambda2[i])})
      
      v1<-sapply(seq(1:M),function(i){var(x1[,i])})
      v2<-sapply(seq(1:M),function(i){var(x2[,i])})
      
      test.stats<-(colSums(x2)/dim(x1)[1]-colSums(x1)/dim(x1)[1])/sqrt(v1/n/Q1+v2/n/(1-Q1))
      e.df<-round((v1/(n*Q1)+v2/(n*(1-Q1)))^2/(v1^2/(n*Q1)^2/(n*Q1-1)+v2^2/(n*(1-Q1))^2/(n*(1-Q1)-1)),0)
      
    }
    t.critical<-qt(1-0.05/2, df=e.df)
    #assurance
    a<-sum(test.stats>t.critical)/M
    
    ### Prior probability that the new treatment is indeed superior
    prior.p<-sum(colSums(x2)>colSums(x1))/M
    
    list(a=a,prior.p=prior.p)
    
  }
  
  est.a<-sapply(N,function(n){assurance(n,x1,x2,v1,v2)})
  a<-unlist(est.a[1,])
  prior.p<-unlist(est.a[2,])
  
  ### Plot power vs assurance
  #dev.new()
  #plot(N,p,type="l",ylim=c(0,1),ylab="Power and assurance",xlab="Total sample size")
  #lines(N,a,type="l",lty=2)
  #legend('bottomright',0.2,c("power", "assurance"),lty=c(1,2))
  
  ### Scatterplots of S_i(1) and S_i(2) for i=1,2 using 100 random samples
  #sampled.surv2<-t(sampled.surv[,1:100])
  #colnames(sampled.surv2)<-c(expression(S1(1)),expression(S1(2)),expression(S2(1)),expression(S2(2)))
  #dev.new()
  #pairs(sampled.surv2)
  
  list(power=p,assurance=a,sample.size=N,prior.p=prior.p)
  
  
}
#example5.2<-assurance.weibull(t1=1,t2=2,N=seq(10,100,20),s11=0.2,s12=0.1,s21=0.3,s22=0.2,
#                              Q1=0.5,M=5000,a.s11=8.10,b.s11=32.81,a.d11=4.36,b.d11=32.61,
#                              m.d12=0.097,v.d12=0.067^2,a.d22=2.23,b.d22=20.10)

inv.ASDWbl <- function(t1,t2,assurance,s11,s12,s21,s22,Q1,M,a.s11,b.s11,a.d11,b.d11,
                       m.d12,v.d12,a.d22,b.d22){
  ##### Functions for estimating the model parameters #####
  
  t1=1;t2=2;assurance=c(0.1,0.2,0.3);s11=0.2;s12=0.1;s21=0.3;s22=0.2;
  Q1=0.5;M=5000;a.s11=8.10;b.s11=32.81;a.d11=4.36;b.d11=32.61;
  m.d12=0.097;v.d12=0.067^2;a.d22=2.23;b.d22=20.10
  e.gamma<-function(s1,s2,t1,t2){
    log(log(s2)/log(s1))/log(t2/t1)
  }
  
  e.lambda<-function(s1,s2,t1,t2){
    -log(s1)/t1^(e.gamma(s1,s2,t1,t2))
  }
  
  e.mu<-function(s1,s2,t1,t2){
    e.lambda(s1,s2,t1,t2)^(-1/e.gamma(s1,s2,t1,t2))*gamma(1+1/e.gamma(s1,s2,t1,t2))
  }
  
  e.sigma2<-function(s1,s2,t1,t2){
    e.lambda(s1,s2,t1,t2)^(-2/e.gamma(s1,s2,t1,t2))*(gamma(1+2/e.gamma(s1,s2,t1,t2))
                                                     -(gamma(1+1/e.gamma(s1,s2,t1,t2)))^2)
  }
  
  
  ### Power function
  powerW2<-function(s11,s12,s21,s22,t1,t2,N,Q1){
    e.mu1<-e.mu(s11,s12,t1,t2)
    e.mu2<-e.mu(s21,s22,t1,t2)
    e.sigma21<-e.sigma2(s11,s12,t1,t2)
    e.sigma22<-e.sigma2(s21,s22,t1,t2)
    v<-e.sigma21/(N*Q1)+e.sigma22/(N*(1-Q1))
    
    powerW2<-pnorm((e.mu2-e.mu1)/sqrt(v)-qnorm(1-0.05/2))
    powerW2
  }
  
  ### Sample s11,s12,s21,s22 for assurance
  
  sampled.s<-function(M,a.s11,b.s11,a.d11,b.d11,m.d12,v.d12,a.d22,b.d22,t1,t2){
    Smatrix<-matrix(NA,nrow=4,ncol=M)
    for(i in 1:M){
      s12<--1;s21<--1;s22<--1;g1<-10;g2<-10
      while(s12<0||s12>1||s21<0||s21>1||s22<0||s22>1||g1>170||g2>170){
        d11<-rbeta(1,a.d11,b.d11)
        d12<-rnorm(1,m.d12,sqrt(v.d12))
        d22<-rbeta(1,a.d22,b.d22)
        s11<-rbeta(1,a.s11,b.s11)
        
        s12<-s11-d11
        s21<-s11+d12
        s22<-s21-d22
        if(s11<0||s12<0){
          g1<-10
        }else{
          g1<-2/e.gamma(s11,s12,t1,t2)
        }
        
        if(s21<0||s22<0){
          g2<-10
        }else{
          g2<-2/e.gamma(s21,s22,t1,t2)
        }
      }
      
      Smatrix[1,i]<-s11
      Smatrix[2,i]<-s12
      Smatrix[3,i]<-s21
      Smatrix[4,i]<-s22
    }
    Smatrix
  }
  
  
  
  ### Calculate power and assuance
  #p<-powerW2(s11,s12,s21,s22,t1,t2,N,Q1)
  sampled.surv<-sampled.s(M,a.s11,b.s11,a.d11,b.d11,m.d12,v.d12,a.d22,b.d22,t1,t2)
  
  ### Simulate the Weibull parameters for two groups
  gamma1<-e.gamma(sampled.surv[1,],sampled.surv[2,],t1,t2)
  lambda1<-e.lambda(sampled.surv[1,],sampled.surv[2,],t1,t2)
  gamma2<-e.gamma(sampled.surv[3,],sampled.surv[4,],t1,t2)
  lambda2<-e.lambda(sampled.surv[3,],sampled.surv[4,],t1,t2)
  
  ### Simulate trial data for each group for each simulated pair of Weibull parameters
  u <- function(n){
      x1<-sapply(seq(1:M),function(i){rweibull(n*Q1,shape=gamma1[i],scale=1/lambda1[i])})
      x2<-sapply(seq(1:M),function(i){rweibull(n*(1-Q1),shape=gamma2[i],scale=1/lambda2[i])})
      
      v1<-sapply(seq(1:M),function(i){var(x1[,i])})
      v2<-sapply(seq(1:M),function(i){var(x2[,i])})
      
      test.stats<-(colSums(x2)/dim(x1)[1]-colSums(x1)/dim(x1)[1])/sqrt(v1/n/Q1+v2/n/(1-Q1))
      e.df<-(v1/(n*Q1)+v2/(n*(1-Q1)))^2/(v1^2/(n*Q1)^2/(n*Q1-1)+v2^2/(n*(1-Q1))^2/(n*(1-Q1)-1))
      
      
      while(sum(e.df==0)>0||sum(is.na(e.df))>0){
        x1<-sapply(seq(1:M),function(i){rweibull(n*Q1,shape=gamma1[i],scale=1/lambda1[i])})
        x2<-sapply(seq(1:M),function(i){rweibull(n*(1-Q1),shape=gamma2[i],scale=1/lambda2[i])})
        
        v1<-sapply(seq(1:M),function(i){var(x1[,i])})
        v2<-sapply(seq(1:M),function(i){var(x2[,i])})
        
        test.stats<-(colSums(x2)/dim(x1)[1]-colSums(x1)/dim(x1)[1])/sqrt(v1/n/Q1+v2/n/(1-Q1))
        e.df<-round((v1/(n*Q1)+v2/(n*(1-Q1)))^2/(v1^2/(n*Q1)^2/(n*Q1-1)+v2^2/(n*(1-Q1))^2/(n*(1-Q1)-1)),0)
        
      }
      t.critical<-qt(1-0.05/2, df=e.df)
      #assurance
      a<-sum(test.stats>t.critical)/M
      
      ### Prior probability that the new treatment is indeed superior
      #prior.p<-sum(colSums(x2)>colSums(x1))/M
      
      #list(a=a,prior.p=prior.p)
      
    }
    
    #a<-getassurance(n)

  
  a = u(300)
  sqrt.inv = inversed(u, lower=0, upper=30000)
  sample <- sqrt.inv(0.43)   
  #new.sample <- list(t0=t0,a=a,b=b,m=m,v=v,M=M,R=R,TT=TT,Q1=Q1,assurance=assurance, sample=sample)
  #est.a<-sapply(N,function(n){assurance(n,x1,x2,v1,v2)})
  #a<-unlist(est.a[1,])
  #prior.p<-unlist(est.a[2,])
  
}
#example5.2<-inv.ASDWbl(t1=1,t2=2,assurance=c(0.1,0.2,0.3),s11=0.2,s12=0.1,s21=0.3,s22=0.2,
 #                             Q1=0.5,M=5000,a.s11=8.10,b.s11=32.81,a.d11=4.36,b.d11=32.61,
  #                            m.d12=0.097,v.d12=0.067^2,a.d22=2.23,b.d22=20.10)

#无删失的ph模型，检验用logrank
#计算所需的所有阳性观察数
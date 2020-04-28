# --------------------------------------------------------------------------------
# ---------—--------------- sample size of bayes Factor --------------------------
# --------------------------------------------------------------------------------


library(BayesFactor)
library(abtest)
library(survival)
library(spBayesSurv)
library(simsurv)
library(doParallel)  
library(BFDA)
# === generate the simulate of data: under H1/BF>1 or H0/BF<1

sample.t.between <- function(n, ES, options.sample=NULL) {
  x <- rnorm(n, mean=0, sd=1)
  y <- rnorm(n, mean=0, sd=1) - ES
  return(cbind(x, y))
}

sample.abtest <- function(n, ES, options.sample=NULL) {
  
  
  
  #增加p1的初始值
  if(is.null(options.sample[["p1"]])){
    p1 <- 0.5
  }else{
    p1 <- as.numeric(options.sample[["p1"]])
  }
  
  if(is.null(options.sample[["effecttype"]])){
    stop("Please define the type of effect size in the argument options.sample.")
  }
  
  effecttype <- options.sample[["effecttype"]]
  
  findp2 <- function(ES, p1, effecttype){
    if(effecttype == "OR"){
      num <- ES*(p1/(1-p1))
      denom <- 1+ES*(p1/(1-p1))
      p2 <- num/denom
    } else if(effecttype == "RR"){
      p2 <- p1*ES
    } else if(effecttype == "AR"){
      p2 <- p1+ES
    }
    
    return(p2)
  }
  if(effecttype == "OR" | effecttype == "logOR"){
    if(effecttype == "logOR") {ES <- exp(ES)} #replace log odds ratio by odds ratio
    p1 <- p1
    p2 <- findp2(ES = ES, p1 = p1, effecttype = "OR")
  }
  
  if(effecttype == "RR"){
    findp1RR <- function(RR){
      i <- 2
      p1 <- p1
      while(1/i * RR >= 1){
        p1 <- 1/(i+1)
        i <- i+1
      }
      return(p1)
    }
    p1 <- findp1RR(RR = ES)
    p2 <- findp2(ES = ES, p1 = p1, effecttype = "RR")
  }
  
  if(effecttype == "AR"){
    if(ES >= 1) stop("Absolute risk can only take values between 0 and 1")
    findp1AR <- function(AR){ # findp1: Find two values whose sum is smaller than one
      i <- 2
      p1 <- p1
      while(1/i + AR >= 1){
        p1 <- 1/(i+2)
        i <- i+1
      }
      return(p1)
    }
    p1 <- findp1AR(AR = ES)
    p2 <- findp2(ES = ES, p1 = p1, effecttype = "AR")
  }
  
  x <- stats::rbinom(n, size = 1, prob = p1)
  y <- stats::rbinom(n, size = 1, prob = p2)
  
  return(cbind(x, y))
}

# === defult the prior is vague prior

# === calculate the BF
BF.t.between <- function(sample, alternative=NULL, prior=NULL) {
  
  #type of hypothesis
  if (alternative=="greater") {
    alt2 <- "BFplus0"
  } else if (alternative=="two.sided"){
    alt2 <- "BF10"
  } else if (alternative=="less"){
    alt2 <- "BFmin0"
  }
  
  
  freq.test <- t.test(sample[, 1], sample[, 2], var.equal=TRUE, alternative=alternative)
  suppressMessages(
    #note the nullInterval:H0：ES=0 vs H1：ES>0
    bf1 <- tryCatch({BayesFactor::ttest.tstat(t = freq.test$statistic,
                                              n1 = nrow(sample), n2 = nrow(sample),
                                              nullInterval=c(0, Inf),
                                              simple=TRUE)},
                    
                    error = function(cond){return(NA)})
    
  )
  #print(prior)
  
  if(prior[[1]] == "Cauchy"){
    
    t1 <- bf10_t(t = as.numeric(freq.test$statistic), n1 = nrow(sample), n2 = nrow(sample), independentSamples = TRUE, prior.location = prior[[2]][["prior.location"]],
                 prior.scale = prior[[2]][["prior.scale"]], prior.df = 1)
    
  } else if (prior[[1]] == "t") {
    
    t1 <- bf10_t(t = as.numeric(freq.test$statistic), n1 = nrow(sample), n2 = nrow(sample), independentSamples = TRUE, prior.location = prior[[2]][["prior.location"]],
                 prior.scale = prior[[2]][["prior.scale"]], prior.df = prior[[2]][["prior.df"]])
    
  } else if (prior[[1]] == "normal") {
    
    
    t1 <- bf10_normal(t = as.numeric(freq.test$statistic), n1=nrow(sample), n2 = nrow(sample), independentSamples = TRUE,
                      prior.mean = prior[[2]][["prior.location"]], prior.variance = prior[[2]][["prior.scale"]])
  }
  
  return(as.numeric(log(t1[[alt2]])))
  # returns the log(BF10)
  #t1 <- bf1[[1]]
  #return(log(t1))
}

BF.abtest <- function(sample, alternative=NULL, freq.test=NULL, prior=NULL, ...) {
  

  if(alternative == "greater") {
    alt2 <- "bfplus0"
  } else if (alternative == "two.sided"){
    alt2 <- "bf10"
  } else if (alternative == "less"){
    alt2 <- "bfminus0"
  }
  
  testdat <- list(y1 = colSums(sample)[1],
                  y2 = colSums(sample)[2],
                  n1 = nrow(sample),
                  n2 = nrow(sample))
  
  testprior <- list(mu_psi = prior[[2]][["prior.mean"]],
                    sigma_psi = prior[[2]][["prior.variance"]],
                    mu_beta = 0,
                    sigma_beta = 1)

  t1 <- as.numeric(abtest::ab_test(data = testdat, prior_par = testprior)$bf[[alt2]])
  #t1 <- as.numeric(abtest::ab_test(data = testdat)$bf[[alt2]])
  return(as.numeric(log(t1)))
  
}
# find the optioal sample size:n
findSSD <- function(Nmin=10, Nmax=1000, power,boundary, ES, alternative,prior=NULL,Tsim=1000, httype, options.sample=NULL){
  
 
  #BF.min <- getn(Nmin)
  #BF.max <- getn(Nmax, ES, alternative)

  Nmax_power <- bfSSD(Nmax, boundary, ES, alternative,prior=prior,Tsim=1000, httype=httype, options.sample=options.sample)$power
  Nmin_power <- bfSSD(Nmin, boundary, ES, alternative,prior=prior,Tsim=1000, httype=httype, options.sample=options.sample)$power
  
  if(Nmax_power<power){
    stop("Nmax is too small")
  }
  
  if(Nmin_power>power){
    stop("Nmin is enough")
  }
  registerDoParallel(cores=4)	
  y<-0
  while(y==0){
    x <- ceiling((Nmin+Nmax)/2)
    x_power <- bfSSD(x, boundary, ES, alternative,prior=prior,Tsim=1000, httype=httype, options.sample=options.sample)$power
    print(c(x,x_power))
    if(Nmax-Nmin==1){
      y <- 1
    }else{
      if(x_power<power ){
        Nmin <- x
      }else{
        Nmax <- x
      }
    }
    
  }
  
  
  alpha <- bfSSD(x, boundary, ES=1, alternative,prior=prior,Tsim=1000, httype=httype, options.sample=options.sample)$power
  

  return(list(sample=x, power=x_power,alpha=alpha))
  
}

bfSSD <- function(N, boundary, ES, alternative,prior=NULL,Tsim=1000, httype, options.sample=NULL){
  
  #Nmin=10;Nmax=100;boundary=3;ES=2;alternative="greater";prior=NULL;Tsim=1000;httype="t"
  if(httype == "t"){
    s.func <- sample.t.between
    bf.func <- BF.t.between
  }else if(httype == "proportion"){
    s.func <- sample.abtest
    bf.func <- BF.abtest
  }
  #note：param change
  getn <- function(n){
    s <- s.func(n, ES, options.sample)
    BF.s <- bf.func(sample=s, alternative=alternative,prior=prior)
    if (BF.s > log(boundary)){
      return(n)
    }else{
      return(0)
    }
  }
  
  
  #calculate the power
  counter <- 0
  for(i in 1:Tsim){
    if(getn(N)>0){
      counter <- counter+1 
    }
  }
  power = counter/Tsim
  
  return(list(sample=N, power=power))
}



#f = findSSD(Nmin=10, Nmax=500, boundary=3, ES=0.12, alternative="greater",prior=NULL,Tsim=1000, httype="proportion", options.sample = options.sample)


BFt <- function(boundary, ES, alternative="greater",priortype="normal", prior_location=0, prior_scale=1, prior_df=1, N,Textparam=list(alternative=c("two.sided", "less", "greater"),priortype=c("normal","t","Cauchy")),CNname=c("阈值","效应量","检验类型","先验类型","位置参数","尺度参数","自由度参数","样本量"),
                tips=list(
                          boundary="预先指定的贝叶斯因子，即所需的证据强度",
                          es= "两组之间的疗效差值",
                          alternative="假设检验的类型，包括双侧(two-sided)，效应量大于0(greater),小于0(less)",
                          priortype="提供关于效应量的三种先验分布柯西(Cauchy),t分布,正态分布",
                          prior_location="先验分布的位置参数",
                          prior_scale="先验分布的尺度参数",
                          prior_df="先验分布的自由度参数(t分布)",
                          N="样本量"
                ),
                defaultValue=list(boundary=6, ES=0.5, alternative="two.sided",priortype="normal", prior_location=0, prior_scale=1, prior_df=0)
                ){
  
  prior = list(priortype,list("prior.location" = prior_location,"prior.scale"= prior_scale, "prior.df"=prior_df))
  Tsim=1000
  #default value
  
  httype="t"
  #prior = NULL
  options.sample = NULL
  
  uu <- function(n){
   getPower <- bfSSD(n, boundary, ES, alternative, prior=prior,Tsim, httype, options.sample)$power
   getalpha <- bfSSD(n, boundary, ES=0, alternative, prior=prior,Tsim, httype, options.sample)$power
   return(c(getPower,getalpha))
  }
  
  powerlist <- mapply(uu,N)
  power <-powerlist[seq(1,length(N)*2-1,2)]
  alpha <-powerlist[seq(2,length(N)*2,2)]
  #sqrt.inv = inversed(u, lower=0, upper=10000)
  #sample <- sqrt.inv(power)   
  conclusion <- paste0("当两组间均数差异为",ES,"，且服从先验分布",priortype,"，分布的位置参数为", prior_location,"，尺度参数为", prior_scale,"，自由度参数为", prior_df,"，贝叶斯因子的阈值为",boundary,"，样本量为",N[1],"时，该条件下的一类错误为",alpha[1],"，功效为",power[1])
  new.sample <- list(boundary=boundary, ES=ES, alternative=alternative,priortype=priortype, prior_location=prior_location, prior_scale=prior_scale, prior_df=prior_df,sample = N, power=power,alpha=alpha,conclusion=conclusion)
  
}



BFpp <- function(N,boundary, ES, alternative="greater",prior_mean=NULL, prior_sd=NULL, effecttype="OR",p1,Textparam=list(alternative=c("two.sided", "less", "greater"),effecttype=c("OR","RR","AR")),CNname=c("样本量","阈值","效应量","检验类型","先验均数","先验方差","效应量类型","发生率"),
                 tips=list(
                   N="样本量",
                   boundary="预先指定的贝叶斯因子，即所需的证据强度",
                   ES= "两组之间的效应量差异值",
                   alternative="假设检验的类型，包括双侧(two-sided)，效应量大于0(greater),小于0(less)",
                   prior_mean="正态分布的第一个参数，相对风险比的对数值log(OR)服从该分布",
                   prior_sd="正态分布的第二个参数，相对风险比的对数值log(OR)服从该分布",
                   effecttype="效应量类型，包括绝对风险(AR),相对风险（RR）OR值",
                   p1 = "组1事件的发生率"
                 ),
                 defaultValue=list(boundary=6, ES=3, alternative="two.sided",prior_mean=0, prior_sd=1, effecttype="OR",p1=0.2)
                 ){
  
  
  #default value
  Tsim=1000
  prior = list("normal",list("prior.mean" = prior_mean,"prior.variance"= prior_sd))
  
  httype="proportion"
  #prior = NULL
  options.sample = c()
  options.sample[["effecttype"]] = effecttype
  options.sample[["p1"]] = p1
  
  uu <- function(n){
    getPower <- bfSSD(n, boundary, ES, alternative, prior,Tsim, httype, options.sample)$power
    if(effecttype=="AR"){
      es0 = 0
    }else{
      es0 = 1
    }
    getalpha <- bfSSD(n, boundary, ES=es0, alternative, prior=prior,Tsim, httype, options.sample)$power
    return(c(getPower,getalpha))
  }
  
  powerlist <- mapply(uu,N)
  power <-powerlist[seq(1,length(N)*2-1,2)]
  alpha <-powerlist[seq(2,length(N)*2,2)]
  #sqrt.inv = inversed(u, lower=0, upper=10000)
  #sample <- sqrt.inv(power)   
  conclusion <- paste0("当相对风险比的对数值为",ES,"，服从先验分布N(",prior_mean,"，",prior_sd,")，组1事件发生率为",p1,"，检验类型为双侧检验，贝叶斯因子的阈值为",boundary,"，样本量为",N[1],"下的I类错误为",alpha[1],"，功效为",power[1])
  new.sample <- list(boundary=boundary, ES=ES, alternative=alternative,prior_mean=prior_mean, prior_sd=prior_sd,effecttype=effecttype,sample = N, p1=p1,power=power,alpha=alpha,conclusion=conclusion)
  
}

####################33
# ==============================================================================
# These are the functions for t-tests with informed priors
# ==============================================================================
# see also https://arxiv.org/abs/1704.02479 for the formulae

#' @import hypergeo

# helper functions for the computation of the Bayes factor with informed priors

A <- function(t, n, nu, mu.delta, g) {
  
  Re(hypergeo::genhypergeo(U = (nu + 1)/2, L = 1/2,
                           z = mu.delta^2*t^2/
                             (2*(1/n + g)*((1 + n*g)*nu + t^2))))
  
}

B <- function(t, n, nu, mu.delta, g) {
  
  out <- mu.delta*t/sqrt(1/2*(1/n + g)*((1 + n*g)*nu + t^2)) *
    exp(lgamma((nu + 2)/2) - lgamma((nu + 1)/2)) *
    Re(hypergeo::genhypergeo(U = (nu + 2)/2, L = 3/2,
                             z = mu.delta^2*t^2/
                               (2*(1/n + g)*((1 + n*g)*nu + t^2))))
  
  return(out)
  
}


C <- function(delta, t, n, nu) {
  
  Re(hypergeo::genhypergeo(U = (nu + 1)/2, L = 1/2,
                           z = n*t^2*delta^2/(2*(nu + t^2))))
  
}

D <- function(delta, t, n, nu) {
  
  out <- t*delta*sqrt(2*n/(nu + t^2))*
    exp(lgamma((nu + 2)/2) - lgamma((nu + 1)/2))*
    Re(hypergeo::genhypergeo(U = (nu + 2)/2, L = 3/2,
                             z = n*t^2*delta^2/(2*(nu + t^2))))
  
  return(out)
  
}

term_normalprior <- function(t, n, nu, mu.delta, g) {
  
  (1 + n*g)^(-1/2) * exp(-mu.delta^2/(2*(1/n + g))) *
    (1 + t^2/(nu*(1 + n*g)))^(-(nu + 1)/2) *
    (A(t, n, nu, mu.delta, g) + B(t, n, nu, mu.delta, g))
  
}

integrand <- function(g, t, n, nu, mu.delta, r, kappa) {
  
  tmp <- term_normalprior(t = t, n = n, nu = nu, mu.delta = mu.delta, g = g)
  pg_log <- kappa/2*(2*log(r) + log(kappa/2)) - lgamma(kappa/2) -
    (kappa/2 + 1)*log(g) - r^2*kappa/(2*g)
  pg <- exp(pg_log)
  out <- tmp*pg
  
  return(out)
  
}

dtss <- function(delta, mu.delta, r, kappa, log = FALSE) {
  
  out <- - log(r) + lgamma((kappa + 1)/2) - .5*(log(pi) + log(kappa)) -
    lgamma(kappa/2) - (kappa + 1)/2 * log(1 + ((delta - mu.delta)/r)^2/kappa)
  
  if ( ! log)
    out <- exp(out)
  
  return(out)
  
}

posterior_t_tmp <- function(delta, t, n1, n2 = NULL, independentSamples = FALSE,
                            prior.location, prior.scale, prior.df,
                            rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, n1*n2/(n1 + n2), n1)
  nu <- ifelse(independentSamples, n1 + n2 - 2, n1 - 1)
  
  mu.delta <- prior.location
  r <- prior.scale
  kappa <- prior.df
  
  numerator <- exp(-neff/2*delta^2)*(1 + t^2/nu)^(-(nu + 1)/2)*
    (C(delta, t, neff, nu) + D(delta, t, neff, nu))*
    dtss(delta, mu.delta, r, kappa)
  
  denominator <- integrate(integrand, lower = 0, upper = Inf,
                           t = t, n = neff, nu = nu, mu.delta = mu.delta,
                           r = r, kappa = kappa, rel.tol = rel.tol)$value
  
  out <- numerator/denominator
  
  if ( is.na(out))
    out <- 0
  
  return(out)
  
}

posterior_t <- Vectorize(posterior_t_tmp, "delta")

cdf_t <- function(x, t, n1, n2 = NULL, independentSamples = FALSE,
                  prior.location, prior.scale, prior.df) {
  
  area <- integrate(posterior_t, lower = -Inf, upper = x, t = t, n1 = n1, n2 = n2,
                    independentSamples = independentSamples,
                    prior.location = prior.location, prior.scale = prior.scale,
                    prior.df = prior.df)$value
  
  if(area > 1)
    area <- 1
  
  return(area)
  
}


posterior_normal_tmp <- function(delta, t, n1, n2 = NULL,
                                 independentSamples = FALSE, prior.mean,
                                 prior.variance,
                                 rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, n1*n2/(n1 + n2), n1)
  nu <- ifelse(independentSamples, n1 + n2 - 2, n1 - 1)
  
  mu.delta <- prior.mean
  g <- prior.variance
  
  numerator <- exp(-neff/2*delta^2)*(1 + t^2/nu)^(-(nu + 1)/2)*
    (C(delta, t, neff, nu) + D(delta, t, neff, nu))*
    dnorm(delta, mu.delta, sqrt(g))
  
  denominator <- term_normalprior(t = t, n = neff, nu = nu,
                                  mu.delta = mu.delta, g = g)
  
  out <- numerator/denominator
  
  if ( is.na(out))
    out <- 0
  
  return(out)
  
}

posterior_normal <- Vectorize(posterior_normal_tmp, "delta")

cdf_normal <- function(x, t, n1, n2 = NULL, independentSamples = FALSE,
                       prior.mean, prior.variance) {
  
  integrate(posterior_normal, lower = -Inf, upper = x, t = t, n1 = n1, n2 = n2,
            independentSamples = independentSamples,
            prior.mean = prior.mean, prior.variance = prior.variance)$value
  
}

# Function to compute the Bayes factor with t distribution as prior

bf10_t <- function(t, n1, n2 = NULL, independentSamples = FALSE, prior.location,
                   prior.scale, prior.df, rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, n1*n2/(n1 + n2), n1)
  nu <- ifelse(independentSamples, n1 + n2 - 2, n1 - 1)
  
  mu.delta <- prior.location
  r <- prior.scale
  kappa <- prior.df
  numerator <- integrate(integrand, lower = 0, upper = Inf,
                         t = t, n = neff, nu = nu, mu.delta = mu.delta,
                         r = r, kappa = kappa,
                         rel.tol = rel.tol)$value
  denominator <- (1 + t^2/nu)^(-(nu + 1)/2)
  
  BF10 <- numerator/denominator
  priorAreaSmaller0 <- integrate(dtss, lower = -Inf, upper = 0,
                                 mu.delta = prior.location, r = prior.scale,
                                 kappa = prior.df)$value
  postAreaSmaller0 <- cdf_t(x = 0, t = t, n1 = n1, n2 = n2,
                            independentSamples = independentSamples,
                            prior.location = prior.location,
                            prior.scale = prior.scale, prior.df = prior.df)
  BFmin1 <- postAreaSmaller0/priorAreaSmaller0
  BFplus1 <- (1 - postAreaSmaller0)/(1 - priorAreaSmaller0)
  BFmin0 <- BFmin1 * BF10
  BFplus0 <- BFplus1 * BF10
  
  return(list(BF10 = BF10, BFplus0 = BFplus0, BFmin0 = BFmin0))
  
}

# Function to compute the Bayes factor with normal distribution as prior

bf10_normal <- function(t, n1, n2 = NULL, independentSamples = FALSE,
                        prior.mean, prior.variance) {
  
  neff <- ifelse(independentSamples, n1*n2/(n1 + n2), n1)
  nu <- ifelse(independentSamples, n1 + n2 - 2, n1 - 1)
  
  mu.delta <- prior.mean
  g <- prior.variance
  numerator <- term_normalprior(t = t, n = neff, nu  = nu,
                                mu.delta = mu.delta, g = g)
  denominator <- (1 + t^2/nu)^(-(nu + 1)/2)
  
  BF10 <- numerator/denominator
  priorAreaSmaller0 <- pnorm(0, mean = prior.mean, sd = sqrt(prior.variance))
  postAreaSmaller0 <- cdf_normal(x = 0, t = t, n1 = n1, n2 = n2,
                                 independentSamples = independentSamples,
                                 prior.mean = prior.mean,
                                 prior.variance = prior.variance)
  BFmin1 <- postAreaSmaller0/priorAreaSmaller0
  BFplus1 <- (1 - postAreaSmaller0)/(1 - priorAreaSmaller0)
  BFmin0 <- BFmin1 * BF10
  BFplus0 <- BFplus1 * BF10
  
  return(list(BF10 = BF10, BFplus0 = BFplus0, BFmin0 = BFmin0))
  
}

#prior_g = list("normal",list("prior.mean" = 1,"prior.variance"= 0.3))

#a = BFt(100,3,0.5,prior.location = 1,prior.scale = 0.3)
#b = BFpp(100,3,1.4,prior.mean = 1,prior.variance= 0.3)


sample.survival.between <- function(N,beta0,lambdas=0.1,gammas=1.5,maxt=5) {
  # Create a data frame with the subject IDs and treatment covariate
  cov <- data.frame(id = 1:(N*2),
                    trt = rbinom(N*2, 1, 0.5))
  
  # Simulate the event times
  dat <- simsurv(lambdas = lambdas, 
                 gammas = gammas, 
                 betas = c(trt = beta0), 
                 x = cov, 
                 maxt = maxt)
  
  # Merge the simulated event times onto covariate data frame
  dat <- merge(cov, dat)
}

BF.survtest <- function(sample,nburn=1000, nsave=1000,nskip=9,prior=prior,dist="weibull"){
  sample["prior"] = 0
  nburn=nburn; nsave=nsave; nskip=nskip; 
  niter=nburn+nsave 
  mcmc=list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=1000) 
  #prior=list(maxL=15,a0=1,b0=1)
  prior=prior; 
  state=list(cpar=1) 
  fit=SuperSurvRegBayes(formula=Surv(eventtime,status)~trt, data=sample,prior=prior,mcmc=mcmc, state=state,dist=dist)
  return(log(fit$BF["PH"]))
}

bfSSDsur <- function(N, boundary, beta0,prior, Tsim=1000){
  #print(prior)
  #Nmin=10;Nmax=100;boundary=3;ES=2;alternative="greater";prior=NULL;Tsim=1000;httype="t"
  #sample.survival.between <- function(N,lambdas=0.1,gammas=1.5,beta0,maxt=5)
  #prior=list(maxL=15,a0=1,b0=1)
  #BF.survtest <- function(sample,nburn=1000, nsave=1000,nskip=9,prior=prior,dist="weibull"){
  
  #note：param change
  getn <- function(n){
    s <- sample.survival.between(n,beta0)
    BF.s <- BF.survtest(s,nburn=1000, nsave=1000,nskip=9,prior=prior,dist="weibull")
    if (BF.s > log(boundary)){
      return(n)
    }else{
      return(0)
    }
  }
  
  
  #calculate the power
  counter <- 0
  for(i in 1:Tsim){
    if(getn(N)>0){
      counter <- counter+1 
    }
  }
  power = counter/Tsim
  
  return(list(sample=N, power=power))
}


BFsur <- function(N, boundary, beta0,prior.maxL=15,prior.a0=0,prior.b0=1, Tsim=1000, Textparam=NULL,CNname=c("样本量","阈值","β0","maxL","a0","b0","模拟次数")){
  
  prior=list(maxL=prior.maxL,a0=prior.a0,b0=prior.b0)
  
  uu <- function(n){
    getPower <- bfSSDsur(n, boundary, beta0,prior, Tsim=Tsim)
    return(getPower$power)
  }
  
  powerlist = as.numeric(lapply(N, uu))
  #sqrt.inv = inversed(u, lower=0, upper=10000)
  #sample <- sqrt.inv(power)   
  new.sample <- list(boundary=boundary, beta0=beta0, prior=prior,Tsim=Tsim, sample = N, level=powerlist)
  
}

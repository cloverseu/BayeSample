library(TrialSize)


OneSampleMean.Equality.classic <- function(alpha, level, sigma, margin,Textparam=NULL,CNname=c("α","power","σ","差异值")){
  sample = OneSampleMean.Equality(alpha, 1-level, sigma, margin)
  new.sample = list(alpha=alpha, level=level, sigma=sigma, margin=margin,sample=sample)
}

TwoSampleMean.Equality.classic <- function(alpha, level, sigma, k, margin,Textparam=NULL,CNname=c("α","power","σ","分配比例","差异值")){
  sample = TwoSampleMean.Equality(alpha, 1-level, sigma,k,margin)
  new.sample = list(alpha=alpha, level=level, sigma=sigma, k=k, margin=margin,sample=sample)
}

TwoSampleProportion.Equality.classic <- function(alpha, level, p1, p2, k, delta,Textparam=NULL,CNname=c("α","power","p1","p2","分配比例","差异值")){
  sample = TwoSampleProportion.Equality(alpha, 1-level, p1, p2, k, delta)
  new.sample = list(alpha=alpha, level=level, p1=p1, p2=p2, k=k, delta=delta,sample=sample)
}

TwoSampleSurvival.Equality.classic <- function(alpha, level, lam1, lam2, k, ttotal, taccrual, gamma, Textparam=NULL,CNname=c("α","power","λ1","λ2","分配比例","总试验时间","实际时间长度","指数分布参数")){
  TwoSampleSurvival.Equality(alpha, 1-level, lam1, lam2, k, ttotal, taccrual, gamma)
  new.sample = list(alpha=alpha, level=level, lam1=lam1, lam2=lam2, k=k, ttotal=ttotal,taccrual=taccrual,sample=sample)
}



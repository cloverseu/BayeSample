library(SampleSizeMeans)
source("D:\\软件\\大论文\\sampleApp\\BayeSample\\classfication\\classicalSize.R", encoding="utf-8")

#code需要优化
##########################
#one group
#1. the variance/mean unkonw
#prior: variance->Gamma(alpha, beta)  mean->norm(no)
#Test for Equality   diff = 0
##########################
mu.varunknow.seq = function(len, alpha, beta, n0, level=0.95, method="acc", Textparam=list(method=c("acc","alc","woc")),CNname=c("长度","α","β","先验样本","概率","方法"),
                            tips=list(len="均值的后验可信区间的长度",
                                   alpha="精度的先验分布（gamma）中的第一个参数(方差的倒数)",
                                   beta="精度的先验分布（gamma）中的第二个参数(方差的倒数)",
                                   n0="等价于均值的先验样本量",
                                   level="后验可信区间的平均覆盖率( average coverage probability ，如0.95)",
                                   method="后验可信区间积分的三种准则(ACC,ALC,WOC)"
                                 ),
                            defaultValue=list(len=0.2,alpha=1,beta=1,n0=1,level=0.95)
                            ){
  
  if (method=="acc"){
    sample  = mu.acc(len, alpha, beta, n0, level)
  }else if (method=="alc"){
    u <- function(level){
      mu.alc(len, alpha, beta, n0, level)
    }
    sample = mapply(u,level)
  }else{
    u <- function(level){
      mu.modwoc(len, alpha, beta, n0, level)
    }
    sample = mapply(u,level)
  }
  #select the ele>1
  conclusion <- paste0("当后验均值的",level[1]*100,"%可信区间长度为",len,"，精度λ~gamma(",alpha,",",
                       beta,")，先验样本量为",n0,"时，采用贝叶斯区间长度法中的",method,"积分准则，试验所需样本量为",sample[1])
  #sample = func(len, alpha, beta, n0, level)
  #sample_c = ceiling(fre.singleNorm(sigma1,0.05,0.2)
  new.sample = list(len=len, alpha=alpha, beta=beta, n0=n0, level=level, method=method, sample=sample,conclusion=conclusion)
}


##########################
#2.varknown:acc,alc,woc are the same
##########################
mu.varknown.seq = function(len, lambda, n0, level = 0.95, Textparam=NULL,CNname=c("长度","λ","先验样本","概率"),
                           tips=list(len="均值的后验可信区间的长度",
                                     lambda="已知的精度 (方差的倒数)",
                                     beta="精度的先验分布（gamma）中的第二个参数(方差的倒数)",
                                     n0="等价于均值的先验样本量",
                                     level="后验可信区间的平均覆盖率( average coverage probability ，如0.95)"
                           ),
                           defaultValue=list(len=0.2,lambda=1,n0=1,level=0.95)
                           ){
  #len = 0.3
  #lambda = 1
  #n0 = c(10,20)
  sample = c()
  listArray = list(len=len, lambda=lambda, n0=n0, level=level)
  if(TRUE %in% (lengths(listArray)>1)){
    seqName = listArray[[names(listArray[lengths(listArray)>1])]]
  }else{
    seqName = 1
  }
  for(i in 1:length(seqName)){
    if (length(len)>1){
      len1 = seqName[i]
    }else{
      len1 = len
    }
    if (length(lambda)>1){
      lambda1 = seqName[i]
    }else{
      lambda1 = lambda
    }
    if (length(n0)>1){
      n1 = seqName[i]
    }else{
      n1 = n0
    }
    if (length(level)>1){
      level1 = seqName[i]
    }else{
      level1 = level
    }
    sample[i] = mu.varknown(len1, lambda1, n1, level1)
  }
  conclusion <- paste0("当后验均值的",level[1]*100,"%可信区间长度为",len,"，精度λ为",lambda,"，先验样本量为",n0,"时，采用贝叶斯区间长度法中，试验所需样本量为",sample[1])
  new.sample = list(len=len, lambda=lambda, level=level, sample=sample,conclusion=conclusion)
}

##########################
#two groups 
#1. the variance unknow,unequal/mean unkonw
#prior: variance->Gamma(alpha1, beta1)/Gamma(alpha2, beta2)  mean->norm(no1,n02)
#Test for Equality   diff = 0
#返回结果为每一组的样本量，允许每一组的样本量不相等
#方差不齐
#mudiff.acc(len=0.2, alpha1=2, beta1=2, alpha2=3, beta2=3, n01=10, n02=25)
##########################
mudiff.unequalvar.seq = function(len, alpha1, beta1, alpha2, beta2, n01, n02, level = 0.95, method="acc",Textparam=list(method=c("acc","alc","woc")),CNname=c("长度","α1","β1","α2","β2","先验样本1","先验样本2","概率","方法"),
                                 tips=list(len="均值的后验可信区间的长度",
                                           alpha1="组1精度的先验分布（gamma）中的第一个参数(方差的倒数)",
                                           beta1="组1精度的先验分布（gamma）中的第二个参数(方差的倒数)",
                                           alpha2="组2精度的先验分布（gamma）中的第一个参数(方差的倒数)",
                                           beta2="组2精度的先验分布（gamma）中的第二个参数(方差的倒数)",
                                           n01="等价于组1 均值的先验样本量",
                                           n02="等价于组2 均值的先验样本量",
                                           level="后验可信区间的平均覆盖率( average coverage probability ，如0.95)",
                                           method="后验可信区间积分的三种准则(ACC,ALC,WOC)"
                                 ),
                                 defaultValue=list(len=0.6, alpha1=1, beta1=1, alpha2=1, beta2=2, n01=10, n02=10)
                                 ){
  #len = 0.3。
  #lambda = 1
  #n0 = c(10,20)
  equal = TRUE
  m = 10000
  mcs = 3
  if (method=="acc"){
    func = mudiff.acc
  }else if (method=="alc"){
    func = mudiff.alc
  }else{
    func = mudiff.modwoc
  }
  
  sample = c()
  listArray = list(len=len, alpha1=alpha1, beta1=beta1, alpha2=alpha2, beta2=beta2, n01=n01, n02=n02, level = level, equal = equal)
  if(TRUE %in% (lengths(listArray)>1)){
    seqName = listArray[[names(listArray[lengths(listArray)>1])]]
  }else{
    seqName = 1
  }
  for(i in 1:length(seqName)){
    if (length(len)>1){
      len1 = seqName[i]
    }else{
      len1 = len
    }
    if (length(alpha1)>1){
      alpha11 = seqName[i]
    }else{
      alpha11 = alpha1
    }
    if (length(alpha2)>1){
      alpha21 = seqName[i]
    }else{
      alpha21 = alpha2
    }
    if (length(beta1)>1){
      beta11 = seqName[i]
    }else{
      beta11 = beta1
    }
    if (length(beta2)>1){
      beta21 = seqName[i]
    }else{
      beta21 = beta2
    }
    if (length(n01)>1){
      n011 = seqName[i]
    }else{
      n011 = n01
    }
    if (length(n02)>1){
      n021 = seqName[i]
    }else{
      n021 = n02
    }
    if (length(level)>1){
      level1 = seqName[i]
    }else{
      level1 = level
    }
    sample[i] = func(len1, alpha11, beta11, alpha21, beta21, n011, n021, level = level1, equal = equal, m = m, mcs = 3)
  }
  conclusion <- paste0("当后验均值的",level[1]*100,"%可信区间长度为",len,"，两组精度分别服从λ1~gamma(",alpha1,",",
                       beta1,"),λ1~gamma(",alpha2,",",beta2,"),","先验样本量分别为",n01,"，",n02,"时，采用贝叶斯区间长度法中的",method,"积分准则，试验所需样本量为",sample[1])
  #sample = func(len, alpha, beta, n0, level)
  new.sample = list(len=len, alpha1=alpha1, beta1=beta1, alpha2=alpha2, beta2=beta2, n01=n01, n02=n02, level = level, equal = equal, method=method, sample=sample,conclusion=conclusion)
}


##########################
#two groups 
#1. the variance unknow,equal/mean unkonw
#prior: variance->Gamma(alpha2, beta2)  mean->norm(no1,n02)
#Test for Equality   diff = 0
#返回结果为每一组的样本量，允许每一组的样本量不相等
#方差不齐
#mudiff.acc(len=0.2, alpha1=2, beta1=2, alpha2=3, beta2=3, n01=10, n02=25)
##########################
mudiff.equalvar.seq <- function(len, alpha, beta, n01, n02, level=0.95, method="acc",Textparam=list(method=c("acc","alc","woc")), CNname=c("长度","α","β","先验样本1","先验样本2","概率","方法"),
                                tips=list(len="均值的后验可信区间的长度",
                                          alpha="精度的先验分布（gamma）中的第一个参数(方差的倒数)",
                                          beta="精度的先验分布（gamma）中的第二个参数(方差的倒数)",
                                          n01="等价于组1 均值的先验样本量",
                                          n02="等价于组2 均值的先验样本量",
                                          level="后验可信区间的平均覆盖率( average coverage probability ，如0.95)",
                                          method="后验可信区间积分的三种准则(ACC,ALC,WOC)"
                                ),
                                defaultValue=list(len=0.6, alpha=1, beta=1,n01=10, n02=10)){
  
  if (method=="acc"){
    func = mudiff.acc.equalvar
  }else if (method=="alc"){
    func = mudiff.alc.equalvar
  }else{
    func = mudiff.modwoc.equalvar
  }
  #select the ele>1
  sample = c()
  listArray = list(len=len, alpha=alpha, beta=beta,n01=n01, n02=n02, level=level)
  if(TRUE %in% (lengths(listArray)>1)){
    seqName = listArray[[names(listArray[lengths(listArray)>1])]]
  }else{
    seqName = 1
  }
  for(i in 1:length(seqName)){
    if (length(len)>1){
      len1 = seqName[i]
    }else{
      len1 = len
    }
    if (length(alpha)>1){
      alpha1 = seqName[i]
    }else{
      alpha1 = alpha
    }
    if (length(n01)>1){
      n011 = seqName[i]
    }else{
      n011 = n01
    }
    if (length(n02)>1){
      n021 = seqName[i]
    }else{
      n021 = n02
    }
    if (length(level)>1){
      level1 = seqName[i]
    }else{
      level1 = level
    }
    if (length(beta)>1){
      beta1 = seqName[i]
    }else{
      beta1 = beta
    }
    sample[i] = func(len1, alpha1, beta1, n011, n021, level1)
  }
  conclusion <- paste0("当后验均值的",level[1]*100,"%可信区间长度为",len,"，两组精度λ~gamma(",alpha,",",
                       beta,")，先验样本量分别为",n01,"，",n02,"时，采用贝叶斯区间长度法中的",method,"积分准则，试验所需样本量为",sample[1])
  #sample = func(len, alpha, beta, n0, level)
  new.sample = list(len=len, alpha=alpha, n01=n01, n02=n02, level=level, method=method, sample=sample,conclusion=conclusion)
}


##########################
#two groups 
#1. the variance unknow,equal/mean unkonw
#prior: variance->Gamma(alpha2, beta2)  mean->norm(no1,n02)
#Test for Equality   diff = 0
#返回结果为每一组的样本量，允许每一组的样本量不相等
#方差不齐
#mudiff.acc(len=0.2, alpha1=2, beta1=2, alpha2=3, beta2=3, n01=10, n02=25)
##########################
#方差已经知道
#mudiff.varknown mudiff.varknown
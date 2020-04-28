library(SampleSizeProportions)


#固定区间概率
propdiff.seq = function(len, c1, d1, c2, d2, level, method="acc",Textparam=list(method=c("acc","alc","woc")),CNname=c("长度","先验1-c1","先验1-d1","先验2-c2","先验2-d2","概率","方法"),
                        tips=list(len="均值的后验可信区间的长度",
                                  c1="beta分布的第一个参数，组1率服从给beta分布",
                                  d1="beta分布的第二个参数，组1率服从给beta分布",
                                  c2="beta分布的第一个参数，组2率服从给beta分布",
                                  d2="beta分布的第二个参数，组2率服从给beta分布",
                                  level="后验可信区间的平均覆盖率( average coverage probability ，如0.95)",
                                  method="后验可信区间积分的三种准则(ACC,ALC,WOC)"
                        ),
                        defaultValue=list(len=0.2, c1=1, d1=1,c2=1, d2=2)
                        ){
    #len = 0.3
    #lambda = 1
    #n0 = c(10,20)
  
    #获取参数的所有变量
    #paramlist <- formalArgs(propdiff.acc.seq)
    equal = TRUE
    if (method=="acc"){
      func = propdiff.acc
    }else if (method=="alc"){
      func = propdiff.alc
    }else{
      func = propdiff.modwoc
    }
    
    sample = c()
    listArray = list(len=len, c1=c1, d1=d1, c2=c2, d2=d2, level = level, equal = equal)
    if(TRUE %in% (lengths(listArray)>1)){
      seqName = listArray[[names(listArray[lengths(listArray)>1])]]
    }else{
      seqName = 1
    }
    #
    for(i in 1:length(seqName)){
      for (j in (names(listArray))){
        if (length(listArray[[j]])>1){
          print(paste0(j,"1"))
          assign(paste0(j,"1"), seqName[i])
        }else{
          print(paste0(j,"1"))
          assign(paste0(j,"1"), listArray[[j]])
        }
        #sample[(i*2-1):(i*2)] = func(len=len1, c1=c11, d1=d11, c2=c21, d2=d21, level = level1, equal = equal, m = 10000, mcs = 3)
      }
      sample[i] = func(len=len1, c1=c11, d1=d11, c2=c21, d2=d21, level = level1, equal = equal1)
    }
    
    conclusion <- paste0("当后验均值的",level[1]*100,"%可信区间长度为",len,"，两组发生率分别服从p1~beta(",c1,",",
                         d1,"),p2~beta(",c2,",",d2,")时，采用贝叶斯区间长度法中的",method,"积分准则，试验所需样本量为",sample[1])
  
    new.sample = list(len=len, c1=c1, d1=d1, c2=c2, d2=d2, level = level, equal = equal, method=method, sample=sample,conclusion=conclusion)
}
  


#mbl方法
propmbl.seq = function(len, c1, d1, c2, d2, level = 0.95, method="acc",Textparam=list(method=c("acc","alc","woc")),CNname=c("长度","先验1-c1","先验1-d1","先验2-c2","先验2-d2","水平","方法")){
  
  if (method=="acc"){
    func = propdiff.mblacc
  }else if (method=="alc"){
    func = propdiff.mblalc
  }else{
    func = propdiff.mblmodwoc
  }
  
  sample = c()
  listArray = list(len=len, c1=c1, d1=d1, c2=c2, d2=d2, level = level)
  if(TRUE %in% (lengths(listArray)>1)){
    seqName = listArray[[names(listArray[lengths(listArray)>1])]]
  }else{
    seqName = 1
  }
  #
  for(i in 1:length(seqName)){
    for (j in (names(listArray))){
      if (length(listArray[[j]])>1){
        print(paste0(j,"1"))
        assign(paste0(j,"1"), seqName[i])
      }else{
        print(paste0(j,"1"))
        assign(paste0(j,"1"), listArray[[j]])
      }
      #sample[(i*2-1):(i*2)] = func(len=len1, c1=c11, d1=d11, c2=c21, d2=d21, level = level1, equal = equal, m = 10000, mcs = 3)
    }
    sample[(i*2-1):(i*2)] = func(len=len1, c1=c11, d1=d11, c2=c21, d2=d21, level = level1, m = 10000, mcs = 3)
  }
  new.sample = list(len=len, c1=c1, d1=d1, c2=c2, d2=d2, level = level,  method=method, sample=sample)
}


#a = propdiff.acc.seq(len=0.05, c1=c(3,4), d1=11, c2=11, d2=54, level = 0.95, equal = FALSE)

library(extraDistr)
library(TailRank)
getNmean <- function(m,y,n1, a,  b){
  succcess <-  a+y
  fail <-  b+n1-y
  
  #calculate the prob of [0,m]
  #ybar <- 0
  #for (i in 0:m){
  #加权平均
  # ybar <- ybar + dbbinom(i, m, alpha = succcess, beta = fail)*i
  #}
  #采用抽样
  ybar <- mean(rbb(100000, m, succcess, fail))
  
  return(ybar)
}
ybarsingle <- function(Tmin,m,y,n1, a,  b){
  print(c(Tmin,m,y,n1, a,  b))
  succcess <-  a+y
  fail <-  b+n1-y
  
  #calculate the prob of [0,m]
  #ybar <- 0
  #for (i in 0:m){
  #加权平均
  # ybar <- ybar + dbbinom(i, m, alpha = succcess, beta = fail)*i
  #}
  #采用抽样
  ybar <- rbbinom(Tmin, m, succcess, fail)
  # if(is.na(ybar)){
  #   print("no")
  #   ybar <- rbbinom(1, m, succcess, fail)
  # }else(
  #   print("yes")
  # )
  #while (is.na(ybar)) {
    #ybar <- rbb(25, m, succcess, fail)[1]
  #}
  #ybar <- rbeta(1,succcess, fail)
  #print(ybar)
  return(ybar)
}

calPower <- function(m,y1,n1,y2,n2,delta,eta, a1,  b1,  a2,  b2, Tmin){
  print(c(m,y1,n1,y2,n2,delta,eta, a1,  b1,  a2,  b2, Tmin))
  #ybar1 <- getNmean(m,y1,n1, a1, b1)
  #ybar2 <- getNmean(m,y2,n2, a2, b2)
  #print(c(ybar1,ybar2))
  #p1 <- rbbinom(Tmin,m,y1+ybar1+ a1, n1+m-y1-ybar1+ b1)/m
  #p2 <- rbbinom(Tmin,m,y2+ybar2+ a2, n2+m-y2-ybar2+ b2)/m
  #到底是哪个的循环抽样
  #pp <- sum(p1-p2>delta)/Tmin
  
  #P1bar <- y1+ybar1/(n1+m)
  #p2bar <- y2+ybar2/(n2+m)
  #pbar <- ((n1+m)*p1bar+(n2+m)*p2bar)/(n1+n2+2*m)
  #Tbin <- (P1bar-p2bar-delta)/((pbar*(1-pbar)*(1/(n1+m)+1/(n2+m)))**0.5)
  #m=5;y1=8;n1=25;y2=3;n2=20;delta=0;eta=0.9; a1=1;  b1=1;  a2=1;  b2=1; Tmin=1000
  
  # count = 0
  # for (i in 1:Tmin){
  #   ybar1 <- ybarsingle(m,y1,n1, a1, b1)
  #   ybar2 <- ybarsingle(m,y2,n2, a2, b2)
  #   p1bar <- (y1+ybar1)/(n1+m)
  #   p2bar <- (y2+ybar2)/(n2+m)
  #   pbar <- ((n1+m)*p1bar+(n2+m)*p2bar)/(n1+n2+2*m)
  #   Tbin <- (p1bar-p2bar-delta)/((pbar*(1-pbar)*(1/(n1+m)+1/(n2+m)))**0.5)
  #   print(c(ybar1,ybar2,p1bar,p2bar,pbar,Tbin))
  #   #print(qnorm(0.95))
  #   print(Tbin>qnorm(0.95))
  #   if(Tbin > qnorm(eta)){
  #     count <- count+1
  #   }
  # }
  # power <- count/Tmin
  
  #b = ppSSR(5:20,20,57,15,59,0,0.95,1,1,1,1,1000)
  
  
    #m=200;y1=20;n1=57;y2=15;n2=59;delta=0;eta=0.9; a1=1;  b1=1;  a2=1;  b2=1; Tmin=10000

    ybar1 <- ybarsingle(Tmin,m,y1,n1, a1, b1)
    ybar2 <- ybarsingle(Tmin,m,y2,n2, a2, b2)
    p1bar <- (y1+ybar1)/(n1+m)
    p2bar <- (y2+ybar2)/(n2+m)
    pbar <- ((n1+m)*p1bar+(n2+m)*p2bar)/(n1+n2+2*m)
    Tbin <- (p1bar-p2bar-delta)/((pbar*(1-pbar)*(1/(n1+m)+1/(n2+m)))**0.5)
    power <- length(Tbin[Tbin>qnorm(0.95)])/length(Tbin)
    
  
  
  #  if(power>eta){
  #   yn <- 1
  #  }else{
  #   yn<-0
  #  }
  # 
  # result <- data.frame(power=power,yn=yn,m=m)
  
}

#未设置所需的最小样本量
ppSSR <- function(N,y1,n1,y2,n2,delta,eta,a1=1, b1=1, a2=1, b2=1,
                  Textparam=NULL,CNname=c("样本量","当前组1成功人数","当前组1总人数","当前组2成功人数","当前组2总人数","δ","η"," a1", " b1", " a2", " b2"),
                  tips=list(
                    N="样本量",
                    y1="当前试验中组1试验成功的人数", 
                    n1="当前试验中组1的总人数",
                    y2="当前试验中组2试验成功的人数", 
                    n2="当前试验中组2的总人数",
                    delta="假设检验中两组间的率差",
                    a1="beta分布的第一个参数，组1率服从beta(a1，b1)",
                    b1="beta分布的第二个参数，组1率服从beta(a1，b1)",
                    a2="beta分布的第一个参数，组2率服从beta(a2，b2)",
                    b2="beta分布的第二个参数，组2率服从beta(a2，b2)",
                    eta="指定率差大于δ时所能接受的试验成功的最小概率，用于有效终止试验"
                  ),
                  defaultValue=list(y1=20,n1=57,y2=12,n2=57,delta=0,eta=0.95,a1=1, b1=1, a2=1, b2=1)
                  ){
  
  #yn_maxm <- calPower(m,y1,n1,y2,n2,delta,eta, a1,  b1,  a2,  b2, Tmin=1000)$yn
  Tmin=1000
  
  u <- function(m){
    calPower(m,y1,n1,y2,n2,delta,eta,a1, b1, a2, b2, Tmin=Tmin)
  }
  power <- mapply(u,N)
  
  conclusion <- paste0("当前试验两组出现试验成功率分别为",y1,"/",n1,"和",y2,"/",n2,"，两组试验成功的先验信息分别为p1~beta(",a1,"，",b1,")，p2~ beta(",a2,"，",b2,")，当前试验δ=",delta,"，η=",eta,",在样本量为",N[1],"，功效为",power[1])
  return(list(sample=N, power=power,conclusion=conclusion))
  #if (yn_maxm == 0){
  # stop("the max sample does't require the condition")
  #}else{
  #优化为最小二乘法
  #}
  
}






#a = calPower(5,8,25,3,20,0.2,0.8,1,1,1,1,1000)

#b = ppSSR(5:20,20,57,15,59,0,0.95,1,1,1,1,1000)



#################################
#PSSD
#################################
library( LearnBayes )

#模拟ys(0->n)
# 分析先验为beta(1,1),似然为ys/n，出现结果大于p0 + delta的后验概率
pre.post.beta <- function( a_a, b_a, p0, delta, n )
{
  ys <- 0:n 
  #分析先验为beta(1,1),似然为ys/n
  pn.s <- (1 - pbeta( p0 + delta, a_a+ys, b_a+n-ys ) )#在这个分布中出现结果大于p0 + delta的概率
  dist <- round (cbind ( n, ys, pn.s ), 4 )
  return ( dist )
}
# 在设计先验中要获得大于s/n的概率有多大
pred.dist <- function( n_d, p_d, a_d, b_d, lambda, n, s, pn.s )
{
  if ( n_d == -1 ) {#若为单点先验，则为二项分布
    md.s <- dbinom( s, n, p_d ) }#在s点的函数值
  else {
    ab_d = c( a_d, b_d )
    md.s <- pbetap( ab_d, n, s ) }
  cum.md.s <- 0
  for ( i in n:s ) {
    if ( n_d == -1 ) {
      wk.md.s <- dbinom( i, n, p_d ) }
    else {
      wk.md.s <- pbetap( ab_d, n, i ) }
    cum.md.s <- cum.md.s + wk.md.s }#算power
  dist <- round( cbind( n, s, pn.s, md.s, cum.md.s ), 4 )
  return( dist )
}

# ???͂??K?v?Ȓl?F
# ???͎??O???zBeta(a,b)?̑?1?????Fa_a
# ???͎??O???zBeta(a,b)?̑?2?????Fb_a
# ?f?U?C?????O???z?̃o???c?L???\???p?????[?^?FnD?i???????AnD = ???̂Ƃ??́An_d = -1?Ƃ????j
# ?f?U?C?????O???z?̃??[?h?Fp_d
# ?ڕW?l?Fp0?{delta?ip0?̓q?X?g???J???ΏƂɂ????鐄???l?Adelta?͏??悹?????̈Ӗ???2?ɕ????Ă??邪?A?ʏ???delta=0?Ƃ????ݒ肪?p????????)
# ?m??臒l?ɁFlambda
# ?m??臒l?��Fgamma
# ?W?{?T?C?Y?v?Z?̎n?_?Fstart?i?ʏ???1?Ƃ????j
# ?W?{?T?C?Y?v?Z?̏I?_?Fend?i?ʏ???200???x?Ƃ????j

#寻找最小样本量
calc_pns_mds <- function( a_a=0.2, b_a=0.8, n_d=-1, p_d=0.5, p0=0.2, delta=0, lambda=0.99, gamma=0.8, start=20, end=100 )
{
  a_d <- 0
  b_d <- 0 
  if ( n_d != -1 ) {
    a_d <- n_d * p_d + 1
    b_d <- n_d * ( 1 - p_d ) + 1
  }
  a <- c(0)
  x <- c(0)
  for ( i in start:end ) {
    #通过分析先验获得要满足max( a[,3] ) > lambda的最小s（/i）是多少
    a <- pre.post.beta ( a_a, b_a, p0, delta, i )
    if ( max( a[,3] ) > lambda ) {
      n <- a[a[,3]==min(a[a[,3]>lambda,3]),1]#n
      s <- a[a[,3]==min(a[a[,3]>lambda,3]),2]#要取得满足max( a[,3] ) > lambda的最小s是多少,s/n
      pn.s <- a[a[,2]==s,3]#后验概率
      #
      b <- pred.dist( n_d, p_d, a_d, b_d, lambda, n, s, pn.s )
      x <- rbind( b, x )
    }
  }
  n <- max( x[x[,5]<gamma,1] ) + 1
  if ( max( x[,5] ) < gamma ) {
    print( 'N is too small !' )
    print( max( x[,5] ) )
  } else { 
    y <- data.frame( 
      n = round( n ), 
      u = round( x[x[,1]==n,2] ), 
      pn.s = round( x[x[,1]==n,3], 4 ), 
      md.s = round( x[x[,1]==n,4], 4 ), 
      pr.d = round( x[x[,1]==n,5], 4 ) ) 
    print( '--- result ---' )
    print( y )
  }
  return(n)
}
#a = calc_pns_mds( a_a=1, b_a=1, n_d=-1, p_d=0.45, p0=0.2, delta=0.1, lambda=0.9)

#由N求power
calc_power <- function( a_a=0.2, b_a=0.8, n_d=-1, p_d=0.5, p0=0.2, delta=0, lambda=0.99, gamma=0.8, i)
{
  
  a_d <- 0
  b_d <- 0 
  if ( n_d != -1 ) {
    a_d <- n_d * p_d + 1
    b_d <- n_d * ( 1 - p_d ) + 1
  }
  a <- c(0)
  x <- c(0)
  
  a <- pre.post.beta ( a_a, b_a, p0, delta, i )
  if ( max( a[,3] ) > lambda ) {
      n <- a[a[,3]==min(a[a[,3]>lambda,3]),1]#n
      s <- a[a[,3]==min(a[a[,3]>lambda,3]),2]#要取得满足max( a[,3] ) > lambda的最小s是多少,s/n
      pn.s <- a[a[,2]==s,3]#后验概率
      #
      b <- pred.dist( n_d, p_d, a_d, b_d, lambda, n, s, pn.s )
  }
  return(b[5])
}
#lambda改成level
#PSSD改成N求power，
PSSD <- function(a_a, b_a,n_d, p_d, p0, delta,lambda,N,Textparam=NULL,CNname=c("分析先验α","分析先验β","设计先验nd","设计先验p","发生率p0","δ","λ"),
                 tips=list(
                    a_a="beta分布的第一个参数，用于分析先验beta(α，β)",
                    b_a="beta分布的第一个参数，用于分析先验beta(α，β)",
                    n_d =" 用于控制设计先验的参数，nd=-1表示nd趋于无穷大，此时设计先验为单点分布",
                    p_d="用于控制设计先验的参数，即备择假设下疗效值的先验信息",
                    p0="对照组或历史数据提供的发生率的参考值",
                    delta="可接受的疗效差异",
                    lambda="分析先验下试验成功的概率",
                    N="样本量"
                 ),
                 defaultValue=list(a_a=0.2, b_a=0.8, n_d=-1, p_d=0.5, p0=0.2, delta=0, lambda=0.99, gamma=0.8)
                 ){
  
  u <- function(m){
    calc_power(a_a=a_a, b_a=b_a, n_d=n_d, p_d=p_d, p0=p0, delta=delta, lambda,i=m)
  }
  power <- mapply(u,N)
  
  conclusion <- paste0("假设检验为H0：θ≤",p0,"+",delta," vs H1：θ≥",p_d,"。θ的分析先验为beta(",a_a,"，",b_a,")，设计先验beta(",nd*pd+1,",",nd(1-pd)+1,")，λ=",lambda,"和γ=",gamma,"的条件下计算试验样本量为",N[1],"power为",power[1])
  return(list(sample=N, power = power, conclusion=conclusion))
  
}

  #yn_maxm <- calPower(m,y1,n1,y2,n2,delta,eta, a1,  b1,  a2,  b2, Tmin=1000)$yn
#a = PSSD( a_a=1, b_a=1, n_d=-1, p_d=0.45, p0=0.2, delta=0.1, level=c(0.9,0.8,0.7))
  #PSSD(a_a=a_a, b_a=b_a, n_d=n_d, p_d=p_d, p0=p0, delta=delta, lambda=m))
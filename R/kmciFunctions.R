### this file contains most of the functions 
### for the kmci package

# uvab: take mean and var of Beta, get back beta parameters
#  u,v may be vectors
uvab<-function(u,v){
    a<-u^2 * (1-u)/v - u
    b<-u*(1-u)^2/v -1+u
    a[u==1 & v==0]<-1
    b[u==1 & v==0]<-0

    a[u==0 & v==0]<-0
    b[u==0 & v==0]<-0
    out<-list(a=a,b=b)
    out
}


# extend qbeta to allow a=0 or b=0 (but not both)
# not needed in R versions >= 3.1.1 (since limits of parameters are defined in qbeta, etc)
# 
if (compareVersion(as.character(getRversion()),"3.1.1")<0){
  qqbeta<-function(x,a,b){
    out<-rep(0,length(a))
    out[b==0]<-1
    out[a==0]<-0
    I<-a>0 & b>0
    if (any(I)) out[I]<-qbeta(x,a[I],b[I])
    # added NA for a==0 and b==0
    # this is for thoroughness, it does not come up in 
    # calculations for bpcp
    out[a==0 & b==0]<-NA
    out
  }
} else {
  # at or after R Version 3.1.1
  # note that qbeta(c(.2,.8),0,0) gives c(0,1)
  qqbeta<-function(x,a,b){ qbeta(x,a,b) }
}

# get bpcp MC from a and b Beta product parameters 
betaprodcimc<-function(a,b,nmc=10^4,alpha=.05){
    np<-length(a)
    # create a matrix of beta distributions
    B<-matrix(rbeta(np*nmc,rep(a,each=nmc),
        rep(b,each=nmc)),nmc,np)
    # multiply columns to create beta product random numbers
    for (j in 2:np){
        B[,j]<-B[,j-1]*B[,j]
    }
    quant<-function(x){ quantile(x,probs=c(alpha/2,1-alpha/2)) }
    apply(B,2,quant)   
}

# get bpcp MM from a and b Beta product parameters
# a,b are vectors of parameters
# returns method of moments estimators for 
#   BP(a[1],b[1]), BP(a[1:2],b[1:2]), ...
betaprodcimm<-function(a,b,alpha=.05){
    ## use notation from Fan,1991
    ## first get Method of moments estimator for the 
    ## cumprod of beta distns
    ## then get associated CIs 
    S<-cumprod( a/(a+b) )
    T<-cumprod( (a*(a+1))/((a+b)*(a+b+1)) )
    cc<- ((S-T)*S)/(T-S^2)
    dd<- ((S-T)*(1-S))/(T-S^2)
    ci<-matrix(NA,2,length(a))
    ci[1,]<-qqbeta(alpha/2,cc,dd)
    ci[2,]<-qqbeta(1-alpha/2,cc,dd)
    ci
}
#betaprodcimc(c(30:20),rep(1,11),nmc=10^6)
#betaprodcimm(c(30:20),rep(1,11))
kmgw.calc<-function(time,status,keepCens=TRUE){
    ## calculate K-M estimator
    ## if keepCens=TRUE then 
    ## output y is all observed times (death or censoring)
    N<-length(time)
    tab<-table(time,status)
    ## if all died or all censored then tab is 
    ## not right dimensions
    dstat<-dimnames(tab)$status
    if (length(dstat)==2){
        di<-tab[,2]
        ci<-tab[,1]
    } else if (length(dstat)==1){
        if (dstat=="1"){
            di<-tab[,1]
            ci<-rep(0,length(tab))
        } else if (dstat=="0"){
            ci<-tab[,1]
            di<-rep(0,length(tab))
        }
    } else stop("status should be 0 or 1")
    
    y<-as.numeric(dimnames(tab)[[1]])
    k<-length(y)

    ni<- c(N,N - cumsum(ci)[-k] - cumsum(di)[-k])
    names(ni)<-names(di)
    ## to avoid overflow problems, make ni, di numeric instead
    ## of interger
    ni<-as.numeric(ni)
    di<-as.numeric(di)
    ci<-as.numeric(ci)
    KM<-cumprod( (ni-di)/ni )
    gw<- KM^2 * cumsum( di/(ni*(ni-di)) )
    ## define variance estimator as 0 after data
    gw[KM==0]<-0

    if (keepCens){
        I<-rep(TRUE,length(y))
    } else {
        I<-di>0
    }
    ## output 
    ##   time=vector of times
    ##     ci[i]=number censored at time[i]
    ##     di[i]=number failed at time[i]
    ##     ni[i]=number at risk just before time[i]
    ##     KM[i]=Kaplan-Meier at time[i]
    ##     gw[i]=Greenwood variance at time[i]
    out<-list(time=y[I],ci=ci[I],di=di[I],ni=ni[I],
        KM=KM[I],gw=gw[I])
    out
}


kmcilog<-function(x,alpha=0.05){
   ### log transformation normal approximation to 
   ###  100(1-alpha) percent CI
   Za<-qnorm(1-alpha/2)
   lower<-exp(log(x$KM) - sqrt(x$gw/x$KM^2) * Za)
   lower[x$KM==0]<-0
   upper<-pmin(1,exp(log(x$KM) + sqrt(x$gw/x$KM^2)*Za))
   upper[is.na(upper)]<-1
   out<-list(time=x$time,surv=x$KM,
       lower=lower,upper=upper,conf.level=1-alpha)
   #class(out)<-"kmci"
   out
}

# intChar creates a character vector describing 
# intervals based on whether L and R is included 
# or not
intChar<-function(L,R,Lin=rep(FALSE,length(L)),
    Rin=rep(TRUE,length(L)),digits=NULL){
    ### check that makes a real interval
    if (length(L)!=length(R)){
        stop("length of L and R must be the same")
    }
    if (any(R-L<0)) stop("L must be less than R")
    Lint<-rep("(",length(L))
    Lint[Lin]<-"["
    Rint<-rep(")",length(L))
    Rint[Rin]<-"]"
    Rint[R==Inf]<-")"
    ### default is to round to the nearest number 
    ### of digits to differentiate 
    ### adjacent times
    if (is.null(digits)){
        x<-sort(unique(c(L,R)))
        md<-min(diff(x))
        if (md<1){
            ### if the minimum distance between adjacent 
            ### times is md
            ### then if we round to digits should 
            ### be able to differentiate
            ### Example:md=.01 then log10(md)=-2 so digits=2
            ###         md=.03 then log10(md)= -1.52 so digits=2
            digits<-ceiling(-log10(md))
        } else {
            digits<-0
        } 
    }
    Interval<-paste(Lint,round(L,digits),",",
        round(R,digits),Rint,sep="")
    Interval
}
#intChar(c(2,4,65),c(3,5,Inf))

# get marks for right censoring times
getmarks<-function(time,status){
    x<-kmgw.calc(time,status)
    ## fix July 29, 2015: get marks when there is 
    ## failure there also
    #x$time[x$ci>0 & x$di==0]
    x$time[x$ci>0]
}
## save time if already have the results from kmgw.calc, use 
getmarks.x<-function(x){
    ## fix July 29, 2015: get marks when there is 
    ## failure there also
    #x$time[x$ci>0 & x$di==0]
    x$time[x$ci>0]
}

borkowf.calc<-function(x,type="log",alpha=.05){
    ## calculate the Borkowf CIs
    ##    type="log" ... do log transformation
    ##    type="logs"...log transformation with shrinkage
    ##    type="norm" ... usual method, no log transformation
    ##    type="norms"...usual method with shrinkage
    ## x is output from kmgw function
    ## take k values and expand for the 2k+1 intervals
    k<-length(x$time)
    d<-rep(x$di,each=2)
    d[2*(1:k)]<-0
    cens<-rep(x$ci,each=2)
    cens[2*(1:k)-1]<-0
    ## add values for interval (0, t1): d=0,n=N,S=1,gw=0
    d<-c(0,d)
    cens<-c(0,cens)
    n<-c(x$ni[1],rep(x$ni,each=2))
    S<-c(1,rep(x$KM,each=2))
    gw<-c(0,rep(x$gw,each=2))

    y<-x$time
    L<-c(0,rep(y,each=2))
    R<-c(rep(y,each=2),Inf)
    Lin<-Rin<-c(rep(c(FALSE,TRUE),k),FALSE)

    #Interval<-intChar(L,R,Lin,Rin)


    #Lint<-rep("(",2*k+1)
    #Lint[Lin]<-"["
    #Rint<-rep(")",2*k+1)
    #Rint[Rin]<-"]"
    #Interval<-paste(Lint,L,",",R,Rint,sep="")

    me<-cumsum(d)
    mc<-cumsum(cens)
    
    n<-x$ni[1]
    bmax<- 1 - me/n
    w<-rep(NA,length(d))

    KM<-S
    w[.5<=KM]<- KM[.5<=KM]
    w[KM<.5 & .5<=bmax]<- .5
    w[bmax<.5]<-bmax[bmax<.5]

    
    VarH<- w*(1-w)/(n - mc)
 
    ## shrunken KM
    KMs<- KM*(1- 1/n) + 1/(2*n)
    ws<- w*(1- 1/n) + 1/(2*n) 
    VarHs<- ws*(1-ws)/(n - mc)

    Za<-qnorm(1-alpha/2)

    lowerNorm<- KM - Za*sqrt(VarH)
    upperNorm<- KM + Za*sqrt(VarH)
    lowerLog<- KM * exp(- Za*sqrt(VarH)/KM )
    upperLog<- KM * exp(+ Za*sqrt(VarH)/KM )

    lowerNorms<- KMs - Za*sqrt(VarHs)
    upperNorms<- KMs + Za*sqrt(VarHs)
    lowerLogs<- KMs * exp(- Za*sqrt(VarHs)/KMs )
    upperLogs<- KMs * exp(+ Za*sqrt(VarHs)/KMs )

    SEG<-sqrt(gw)
    SEH<-sqrt(VarH)
    SEHs<-sqrt(VarHs)

    n.minus.mc<-n-mc

    fixup<-function(x){ x[x>1 | is.na(x)]<-1; x }
    fixlo<-function(x){ x[x<0 | is.na(x)]<-0;x }

    lowerNorm<-fixlo(lowerNorm)
    lowerNorms<-fixlo(lowerNorms)
    lowerLog<- fixlo(lowerLog)
    lowerLogs<-fixlo(lowerLogs)
 
    upperNorm<-fixup(upperNorm)
    upperNorms<-fixup(upperNorms)
    upperLog<-fixup(upperLog)
    upperLogs<-fixup(upperLogs)

  

    #out<-data.frame(time,n.minus.mc,SEG,SEH,SEHs,
    #     ci,di,ni,me,mc,KM,bmax,w,KMs,ws,VarH,VarHs,
    #     lowerNorm,upperNorm,lowerLog,upperLog,
    #     lowerNorms,upperNorms,lowerLogs,upperLogs)


    #    out<-list(y=x$time,d=d,n=n,S=S,gw=gw)
    if (type=="log"){
         surv<-KM
        lower<-lowerLog
        upper<-upperLog
    } else if (type=="logs"){
        surv<-KMs
        lower<-lowerLogs
        upper<-upperLogs
    } else if (type=="norm"){
        surv<-KM
        lower<-lowerNorm
        upper<-upperNorm
    } else if (type=="norms"){
        surv<-KMs
        lower<-lowerNorms
        upper<-upperNorms
    }

    #out<-data.frame(L=L,Lin=Lin,R=R,Rin=Rin,
    #    SEG=SEG,SEH=SEH,SEHs=SEHs,
    #    d=d,cens=cens,n=n,surv=surv,KM=KM,
    #    me=me,mc=mc,bmax=bmax,w=w,KMs=KMs,
    #    ws=ws,VarH=VarH,VarHs=Vars,gw=gw,
    #    lowerNorm=lowerNorm,upperNorm=upperNorm,
    #    lowerLog=lowerLog,upperLog=upperLog,
    #    lowerNorms=lowerNorms,upperNorms=upperNorms,
    #    lowerLogs=lowerLogs,
    #    upperLogs=upperLogs,lower=lower,
    #    upper=upper,row.names=Interval)
    #
    # getmarks.x gets censoring marks for plotting
    out<-list(cens=getmarks.x(x),L=L,Lin=Lin,R=R,Rin=Rin,
        surv=surv,lower=lower,upper=upper,conf.level=1-alpha)
    out
}

kmciBorkowf<-function(time,status,type="log",alpha=0.05){
    x<-kmgw.calc(time,status,keepCens=TRUE)
    out<-borkowf.calc(x,type,alpha)
    class(out)<-"kmciLR"
    out
}

kmci1TG<-function(time,status,tstar,alpha=.05){
    x<-kmgw.calc(time,status,keepCens=FALSE)

    I<- x$time <= tstar & (x$ni>x$di)
    if (any(I)){

    ni<-x$ni[I]
    di<-x$di[I]
    R<-function(lambda){
         sum( (ni-di)*log(1+lambda/(ni-di)) - 
                ni*log(1+lambda/ni) )
    }

    Chisq<- qchisq(1-alpha,1)
    rootfunc<-function(lambda){
       nlam<-length(lambda)
       out<-rep(NA,nlam)
       for (i in 1:nlam){
           out[i]<- -2*R(lambda[i]) - Chisq
       }
       out
    }

    lam1<-uniroot(rootfunc,c(0,10^6))$root
    lam2<-uniroot(rootfunc,c(-min(ni-di),0))$root
    Supper<- prod( (ni+lam1-di)/(ni+lam1) )
    Slower<- prod( (ni+lam2-di)/(ni+lam2) )
    KM<-prod( (ni-di)/(ni) )
    } else {
    ### no time<tstar
        KM<-1
        Supper<-1
        Slower<-1
    }
    out<-list(time=tstar,surv=KM,upper=Supper,lower=Slower)
    class(out)<-"kmci"
    out
}

kmciTG<-function(time,status,alpha=.05){
    utime<-sort(unique(time[status==1]))
    ntime<-length(utime)
    surv<-upper<-lower<-rep(NA,ntime)
    for (i in 1:ntime){
        x<-kmci1TG(time,status,utime[i],alpha)
        surv[i]<-x$surv
        lower[i]<-x$lower
        upper[i]<-x$upper
    }
    ## if largest observation is censored, 
    ## then need to add on extra value at the end
    if (max(time)>max(utime)){
        tmax<-max(time)
        k<-length(surv)
        out<-list(cens=getmarks(time,status),
          time=c(utime,tmax),surv=c(surv,surv[k]),
          upper=c(upper,upper[k]),
          lower=c(lower,lower[k]),conf.level=1-alpha)
    } else {
        out<-list(cens=getmarks(time,status),
          time=utime,surv=surv,upper=upper,
          lower=lower,conf.level=1-alpha)
    }
    class(out)<-"kmci"
    out
}


kmConstrain<-function(tstar,pstar,x,alpha=.05){
    ## KM given that S(t)=pstar
    ## first find index such that
    ## 
    I<- x$time <= tstar
    if (any(I)){
        nj<-x$ni[I]
        dj<-x$di[I]
        rootfunc<-function(lambda){
            nl<-length(lambda)
            out<-rep(NA,nl)
            for (i in 1:nl){
                out[i]<- prod( (nj + lambda[i] - dj)/
                                (nj+lambda[i]) ) - pstar
            }
            out
        }
        lambda<-uniroot(rootfunc,c(-min(nj-dj),
                              10^3 * nj[1]))$root    
        ## now calculate constrained variance
        pbar<-(nj + lambda -dj)/(nj+lambda)
        Sbar<- cumprod( pbar )
        out<-list(time=x$time[I],Sc=Sbar)
    } else {
        out<-list(time=tstar,Sc=pstar)
    }
    out
}
rejectFromInt<-function(theta,interval,thetaParm=FALSE){
    ## thetaParm=TRUE means theta is the true value, 
    ## and interval is a confidence interval
    ## thetaParm=FALSE means theta is an estimate and 
    ## interval is the null distribution
    if (length(interval)!=2) stop("interval should be length 2")
    int<-c(min(interval),max(interval))
    reject<-rep(0,3)
    names(reject)<-c("estGTnull","estLTnull","two.sided")
    ## if thetaParm=TRUE then theta is the parameter 
    ## under the null and 
    ## interval is a confidence interval
    ## so if theta<int[1], then 
    ## the estimate is greater than the null hypothesized 
    ## value of the parameter
    if (thetaParm){
        reject["estGTnull"]<- ifelse(theta<int[1],1,0)
        reject["estLTnull"]<- ifelse(theta>int[2],1,0)
    } else {
    ## if thetaParm=FALSE then theta is an estimate and 
    ## interval are quantile from a null distribution
    ## so if theta<int[1], then 
    ## the estimate is less than the null hypothesized 
    ## value of the parameter
        reject["estGTnull"]<- ifelse(theta>int[2],1,0)
        reject["estLTnull"]<- ifelse(theta<int[1],1,0)
    }
    ### two-sided rejection is sum of one sided rejections
    reject[3]<-sum(reject)
    reject
}

kmtestBoot<-function(time,status,tstar,pstar,M=1000,alpha=0.05){
    x<-kmgw.calc(time,status,keepCens=FALSE)
    ### pick out KM at tstar
    Sx<-function(x,tstar){
        I<-x$time<=tstar
        if (!any(I)) out<-1
        else out<-x$KM[max((1:length(x$KM))[I])]
        out
    }

    # get observed value
    Sobs<-Sx(x,tstar)
 
    n<-length(time)

    SB<-rep(NA,M)
    for (i in 1:M){
        ii<-sample(n,replace=TRUE)   
        temp<-kmgw.calc(time[ii],status[ii],keepCens=FALSE)
        SB[i]<-Sx(temp,tstar)
    }
    ### use type=4 quantile so that equals value defined in
    ### Barber and Jennison  S[M*0.025] = S[25] when M=1000
    quantilesNullDistribution<-quantile(SB,c(alpha/2,
                             1-alpha/2),type=4)
    reject<-rejectFromInt(pstar,quantilesNullDistribution,
              thetaParm=TRUE)
    reject
}

kmtestConstrainBoot<-function(time,status,tstar,pstar,M=1000,alpha=0.05){
    x<-kmgw.calc(time,status,keepCens=FALSE)
    ### pick out KM at tstar
    Sx<-function(x,tstar){
        I<-x$time<=tstar
        if (!any(I)) out<-1
        else out<-x$KM[max((1:length(x$KM))[I])]
        out
    }

    Sobs<-Sx(x,tstar)

    xcon<-kmConstrain(tstar,pstar,x)
    ## calculate censoring distribution by reversing status
    xcens<-kmgw.calc(time,1-status,keepCens=FALSE)
    n<-length(time)
    xSurv<-c(xcon$time,Inf)
    ## get density function, add extra element 
    ## for all survival times after tstar
    dSurv<- -diff(c(1,xcon$Sc))
    dSurv<- c(dSurv,1-sum(dSurv))
    ## in case last element is death, so KM of 
    ## censoring distribution does not 
    ## go to zero, add extra element at max(time)+1
    xCens<- c(xcens$time,max(time)+1)
    dCens<- -diff(c(1,xcens$KM))
    dCens<-c(dCens,1-sum(dCens))
    SB<-rep(NA,M)
    for (i in 1:M){
        if (length(dCens)>1){ 
           Ci<-sample(xCens,n,prob=dCens,replace=TRUE)
        } else Ci<-rep(xCens,n)
        Xi<-sample(xSurv,n,prob=dSurv,replace=TRUE)
        Time<-Xi
        Time[Xi>Ci]<-Ci[Xi>Ci]
        StatusTF<- Xi==Time
        Status<-rep(0,n)
        Status[StatusTF]<-1
        temp<-kmgw.calc(Time,Status)
        ### pick out KM at tstar
        SB[i]<-Sx(temp,tstar)
    }
    ### use type=4 quantile so that equals value defined 
    ### in Barber and Jennison  S[M*0.025] = S[25] when M=1000
    quantilesNullDistribution<-quantile(SB,c(alpha/2,
        1-alpha/2),type=4)
    reject<-rejectFromInt(Sobs,quantilesNullDistribution,
        thetaParm=FALSE)
    reject
}


kmConstrainBeta.calc<-function(tstar,pstar,x,alpha=.05){
    ## KM given that S(t)=pstar
    ## first find index such that
    ## 
    if (length(tstar)>1) stop("tstar must be a scalar")
    I<- x$time <= tstar
    nj<-x$ni[I]
    dj<-x$di[I]
    if (all(I==FALSE)){
        ### no x$time before tstar
        ### then dj=0, and it does not matter 
        ### what nj and lamba are since when dj=0
        ### pbar=1 for all lambda
        ### and qbar=0
        ### so vc=0 and just define Sobs=1 (KM before 
        ###     first death), 
        ### lower=1 and upper=1
        Sobs<-1
        qlower<-1
        qupper<-1        
    } else {
        rootfunc<-function(lambda){
            pstar - prod( (nj + lambda - dj)/(nj+lambda) )
        }
        lambda<-uniroot(rootfunc,c(-min(nj-dj),
             10^3 * nj[1]))$root    
        ## now calculate constrained variance
        pbar<-(nj + lambda -dj)/(nj+lambda)
        qbar<-1-pbar
        Sbar<- cumprod( pbar )
        Shat<-x$KM[I]
        ## use Sbar and Shat just before tj, so add 1 
        ## to beginning and delete last
        ns<-length(Sbar)
        Sbar<-c(1,Sbar[-ns])
        Shat<-c(1,Shat[-ns])
        vc<- (pstar^2)*sum(  (Shat*qbar)/(nj*Sbar*pbar) ) 
        abc<-uvab(pstar,vc)
        ### beta confidence interval
        qlower<-qqbeta(alpha/2,abc$a,abc$b)
        qupper<-qqbeta(1-alpha/2,abc$a,abc$b)
        ### pick out KM at tstar
        Sx<-function(x,tstar){
            I<-x$time<=tstar
            if (!any(I)) out<-1
            else out<-x$KM[max((1:length(x$KM))[I])]
            out
         }
         Sobs<-Sx(x,tstar)
    }     
    quantilesNullDistribution<-c(qlower,qupper)
    reject<-rejectFromInt(Sobs,quantilesNullDistribution,
      thetaParm=FALSE)
    reject
}

kmtestConstrainBeta<-function(time,status,tstar,pstar,alpha=.05){
     x<-kmgw.calc(time,status,keepCens=FALSE)
     out<-kmConstrainBeta.calc(tstar,pstar,x,alpha)
     out
}
#x<-kmgw.calc(leuk$time,leuk$status,keepCens=FALSE)
#kmtestConstrainBeta(leuk$time,leuk$status,10,.5)
#kmtestConstrainBoot(leuk$time,leuk$status,10,.5)

kmciSW<-function(time,status,alpha=.05){
    ## This function gives confidence intervals following 
    ## Strawderman and Wells (1997, JASA, 1356-1374)
    ## notation follows
    ## Strawderman, Parzen, and Wells (1997, 
    ## Biometrics 1399-1415
    ## 
    ## for this method we need the Nelson-Aalen estimator
    x<-kmgw.calc(time,status,keepCens=FALSE)
    ## here we allow grouped survival data, and 
    ## use the usual N-A estimator for it
    Lambda<-cumsum(x$di/x$ni)
    ## eq 2: sig= sigmahat_A
    sig<- sqrt(cumsum(1/x$ni^2))
    ## eq 10: kappa
    kappa<- (1/sig^3) * cumsum(1/x$ni^3)
    ## eq 9: we need to use this twice, once for the 
    ## lower interval, once for upper
    ## first lower interval (upper for Lambda, lower for S)
    Za<- qnorm(alpha/2)
    LambdaSWlower<-Lambda - sig*(Za + 
        (sig/4 - kappa/3)*Za^2 - (sig/4+kappa/6))
    Slower<-exp(-LambdaSWlower)
    ## now upper (lower for Lambda, upper for S)
    Za<- qnorm(1-alpha/2)
    LambdaSWupper<-Lambda - sig*(Za + 
        (sig/4 - kappa/3)*Za^2 - (sig/4+kappa/6))
    Supper<-exp(-LambdaSWupper)

    ## if largest observation is censored, then need to add
    ##  on extra value at the end
    ## so that the KM plots all the way until the last 
    ## censored observation
    ##  lower and upper also stay the same out to the 
    ## last censoring value
    if (max(time)>max(x$time)){
        tmax<-max(time)
        k<-length(x$time)
        x$time<-c(x$time,tmax)
        x$KM<-c(x$KM,x$KM[k])
        Slower<-c(Slower,Slower[k])
        Supper<-c(Supper,Supper[k])
    } 
    out<-list(cens=getmarks(time,status),time=x$time,
       surv=x$KM,lower=Slower,upper=Supper,
       conf.level=1-alpha,
       ni=x$ni,di=x$di,Lambda=Lambda,SNA=exp(-Lambda),
       LambdaLower=LambdaSWlower,LambdaUpper=LambdaSWupper)

    class(out)<-"kmci"
    out
}


kmtestBinomial<-function(time,status,cens,t0,S0,alpha=0.05){
    ## test H0: S(t0)=S0
    ## using exact binomial test with only individuals 
    ## who had censoring times after t0
    ##
    I<- cens>t0
    ## this is a slow way to do it, but it is easier to program
    out<-bpcp(time[I],status[I],nmc=0,alpha=alpha,
        Delta=0,stype="km")
    class(out)<-"kmciLR"
    sci<-StCI(out,t0)
    out<-rejectFromInt(S0,sci[3:4],thetaParm=TRUE)
    out
}


kmtestALL<-function(time,status,t0,S0,cens=NULL,M=1000,NMC=10^5,alpha=0.05){

    maxDeath.or.Cens<-max(time)

    ### find reject for normal log transform method
    rnormlog<-function(){
        x<-kmcilog(kmgw.calc(time,status,keepCens=FALSE),alpha)
        sci<-StCI(x,tstar=t0)
        reject<-rejectFromInt(S0,c(sci$lower,sci$upper),
            thetaParm=TRUE)
        reject
    }

    ### Strawderman-Wells method
    rSW<-function(){
        x<-kmciSW(time,status,alpha)
        sci<-StCI(x,t0)
        reject<-rejectFromInt(S0,c(sci$lower,sci$upper),
            thetaParm=TRUE)
        reject
    }

    ### Borkowf method 
    rBlog<-function(){
        x<-kmciBorkowf(time,status,type="log",alpha)
        sci<-StCI(x,t0)
        reject<-rejectFromInt(S0,c(sci$lower,sci$upper),
            thetaParm=TRUE)
        reject
    }

    ### Borkowf method shrinkage 
    rBlogs<-function(){
        x<-kmciBorkowf(time,status,type="logs",alpha)
        sci<-StCI(x,t0)
        reject<-rejectFromInt(S0,c(sci$lower,sci$upper),
            thetaParm=TRUE)
        reject
    }


    ### constrained bootstrap method
    rcboot<-function(){
        kmtestConstrainBoot(time,status,t0,S0,M,alpha=alpha)
    }

    ### constrained beta method 
    rcbeta<-function(){
        kmtestConstrainBeta(time,status,t0,S0,alpha)
    }

    ### binomial method 
    rbinom<-function(){
        kmtestBinomial(time,status,cens,t0,S0,alpha)
    }

    ### likelihood ratio (Thomas and Grunkemeier) method 
    rTG<-function(){
        sci<-kmci1TG(time,status,t0,alpha)
        reject<-rejectFromInt(S0,c(sci$lower,sci$upper),
            thetaParm=TRUE)
        reject
    }

    ### Repeated Beta method with method of moments 
    rRBmm<-function(){
        x<-bpcp(time,status,nmc=0,alpha)
        sci<-StCI(x,t0)
        reject<-rejectFromInt(S0,c(sci$lower,sci$upper),
            thetaParm=TRUE)
        reject
    }


    ### Repeated Beta method with Monte Carlo 
    rRBmc<-function(){
        x<-bpcp(time,status,nmc=NMC,alpha=.05)
        sci<-StCI(x,t0)
        reject<-rejectFromInt(S0,c(sci$lower,sci$upper),
            thetaParm=TRUE)
        reject
    }

    ## get names for rejection values from rnormlog()
    normlog<-rnormlog()
    Rejections<-matrix(NA,10,3,dimnames=list(
        c("normlog","SW","Blog","Blogs","cbeta",
          "TG","cboot","binom","RBmm","RBmc"),
        names(normlog)))
    Rejections["normlog",]<-normlog
    Rejections["SW",]<-rSW()
    Rejections["Blog",]<-rBlog()
    Rejections["Blogs",]<-rBlogs()
    Rejections["cbeta",]<-rcbeta()
    Rejections["TG",]<-rTG()
    Rejections["RBmm",]<-rRBmm()
    if (!is.null(cens)) Rejections["binom",]<-rbinom()
    if (M>0) Rejections["cboot",]<-rcboot()
    if (NMC>0) Rejections["RBmc",]<-rRBmc()
    Rejections
}

#kmtestALL(1:20,rep(1,20),10,.9)

getDefault.mark.time<-function(inmt,inx){
    ### if mark.time=NULL then: stype="mue" then do 
    ### not use mark.time, if stype="km" use mark.time
    if (is.null(inmt)){
       if (is.null(inx$stype)){ inmt<-TRUE
       } else inmt<- inx$stype=="km"
    }
    inmt
}

plot.kmciLR<-function(x,XLAB="time",YLAB="Survival",
    YLIM=c(0,1),ciLTY=2,
    ciCOL=gray(.8),mark.time=NULL,linetype="both",...){

    mark.time<-getDefault.mark.time(mark.time,x)

    ### if xlab, ylab, or ylim is NOT given explicitly, 
    ### then replace with default XLAB, YLAB or YLIM
    ### so we need to get the ... from the call
    ### if xlab, ylab or ylim are not in ... then add 
    ### the default
    ### then when we do the plot, we use the call function 
    ### so we can use anything 
    ### else that was in the ...
    md<-match.call(expand.dots=FALSE)
    dots<-as.list(md)[["..."]]
    class(dots)<-"list" 
    if (is.null(dots[["xlab"]])) dots$xlab<-XLAB
    if (is.null(dots[["ylab"]])) dots$ylab<-YLAB
    if (is.null(dots[["ylim"]])) dots$ylim<-YLIM
    do.call("plot",c(list(
        x=c(x$L,x$R[x$R<Inf]),
        y=c(x$surv,x$surv[x$R<Inf]),
        type="n"), dots) )

    #plot(c(x$L,x$R[x$R<Inf]),c(x$surv,x$surv[x$R<Inf]),
    #    ylim=YLIM,type="n",xlab=XLAB,ylab=YLAB,...)
    n<-length(x$L)

    if (linetype=="both" | linetype=="surv"){
        ## horizontal KM
        segments(x$L,x$surv,x$R,x$surv,...)
        ## vertical lines connecting KM estimates
        segments(x$L[-1],x$surv[-1],x$R[-n],x$surv[-n],...)
        ### add crosses for censored objects
        if (mark.time){
            xcens<-x$cens
            if (is.null(xcens)){
                ## get censoring times when x$cens does 
                ## not exist
                if (!is.null(x$n.censor) & !is.null(x$time)){
                    xcens<- x$time[x$n.censor>0]       
                } else {
                    warning("censoring times not plotted")
                    xcens<-numeric(0)
                }
            } 
            ## if no censoring times, then xcens=numeric(0)
            if (length(xcens)>0){
                out<-StCI(x,xcens)
                points(out$time,out$survival,pch=3)
            }
        }
    }
    if (linetype=="both" | linetype=="ci"){
        segments(x$L,x$lower,x$R,
            x$lower,lty=ciLTY,col=ciCOL,...)
        segments(x$L[-1],x$lower[-1],x$R[-n],
            x$lower[-n],lty=ciLTY,col=ciCOL,...)

        segments(x$L,x$upper,x$R,
            x$upper,lty=ciLTY,col=ciCOL,...)
        segments(x$L[-1],x$upper[-1],x$R[-n],
            x$upper[-n],lty=ciLTY,col=ciCOL,...)
    }

}

plot.kmciLRtidy <- function(x, ...) {
  tidyout <- tidykmciLR(x)
  if (ncol(tidyout) == 5) {
    ggplot(tidyout, aes_string(x = "time", y = "surv", ymin = "lower", ymax = "upper", linetype = "group")) + 
      geom_line() + geom_ribbon(alpha = .2) + labs(linetype = x[[1]]$groupVarName) + xlab("Time") + ylab("Survival")
  }
  else {
    ggplot(tidyout, aes_string(x = "time", y = "surv", ymin = "lower", ymax = "upper")) + 
      geom_line() + geom_ribbon(alpha = .2) + xlab("Time") + ylab("Survival")
  }
}


plot.kmciLRgroup <- function(x,XLAB="Time",YLAB="Survival",
                             YLIM=c(0,1),ciLTY=2,
                             ciCOL=gray(.8), linetype="both",...) {
  tidyout <- tidykmciLR(x)
  md<-match.call(expand.dots=FALSE)
  dots<-as.list(md)[["..."]]
  class(dots)<-"list" 
  if (is.null(dots[["xlab"]])) dots$xlab<-XLAB
  if (is.null(dots[["ylab"]])) dots$ylab<-YLAB
  if (is.null(dots[["ylim"]])) dots$ylim<-YLIM
  do.call("plot",c(list(
    x=c(tidyout$time,tidyout$time[tidyout$time<Inf]),
    y=c(tidyout$surv,tidyout$surv[tidyout$time<Inf]),
    type="n"), dots) )
  
  #plot(c(x$L,x$R[x$R<Inf]),c(x$surv,x$surv[x$R<Inf]),
  #    ylim=YLIM,type="n",xlab=XLAB,ylab=YLAB,...)
  n<-length(tidyout$time)
  
  #If a group variable is provided, change the line width based on the level of the treatment variable
  #and provide a legend
  #If no group variable is provided, don't change the line width
  if (linetype=="both" | linetype=="surv"){
    segments(tidyout$time,tidyout$surv,tidyout$time,tidyout$surv, lwd=as.numeric(factor(tidyout$group)), ...)
    segments(tidyout$time[-1],tidyout$surv[-1],tidyout$time[-n],tidyout$surv[-n], lwd=as.numeric(factor(tidyout$group)), ...)
    legend("topright", legend = unique(tidyout$group), lwd = unique(as.numeric(factor(tidyout$group))))
  }
  
  if (linetype=="both" | linetype=="ci"){
    segments(tidyout$time,tidyout$lower,tidyout$time,
             tidyout$lower,lty=ciLTY,col=ciCOL,...)
    segments(tidyout$time[-1],tidyout$lower[-1],tidyout$time[-n],
             tidyout$lower[-n],lty=ciLTY,col=ciCOL,...)
    
    segments(tidyout$time,tidyout$upper,tidyout$time,
             tidyout$upper,lty=ciLTY,col=ciCOL,...)
    segments(tidyout$time[-1],tidyout$upper[-1],tidyout$time[-n],
             tidyout$upper[-n],lty=ciLTY,col=ciCOL,...)
  }
}


citoLR<-function(x){
    ## add L and R to kmci object
    x$L<-c(0,x$time)
    x$R<-c(x$time,Inf)
    x$Lin<-rep(FALSE,length(x$time)+1)
    x$Rin<-c(rep(TRUE,length(x$time)),FALSE)
    x$surv<-c(1,x$surv)
    x$lower<-c(1,x$lower)
    x$upper<-c(1,x$upper)
    x
}
plot.kmci<-function(x,...){
    xLR<-citoLR(x)
    ### because plot.kmciLR calls StCI, need to make sure call StCI.kmciLR, not default
    class(xLR)<-"kmciLR"
    plot.kmciLR(xLR,...)
}

lines.kmci<-function(x,...){
    xLR<-citoLR(x)
    if (!is.null(xLR$cens)){
        xLR$R[xLR$R==Inf]<-max(xLR$cens)
    } else {
        xLR$R[xLR$R==Inf]<-max(xLR$time)
    }
    lines.kmciLR(xLR,...)
}


lines.kmciLR<-function(x,lty=c(2,1),col=c(gray(.8),gray(0)),linetype="both",mark.time=NULL,...){
    if (length(lty)==1) lty<-rep(lty,2)
    mark.time<-getDefault.mark.time(mark.time,x)
    n<-length(x$L)
    if (linetype=="ci" | linetype=="both"){
        segments(x$L,x$lower,x$R,x$lower,
            lty=lty[1],col=col[1],...)
        segments(x$L[-1],x$lower[-1],x$R[-n],x$lower[-n],
            lty=lty[1],col=col[1],...)
        segments(x$L,x$upper,x$R,x$upper,
            lty=lty[1],col=col[1],...)
        segments(x$L[-1],x$upper[-1],x$R[-n],x$upper[-n],
            lty=lty[1],col=col[1],...)
    }
    if (length(lty)==1) lty<-rep(lty,2)
    if (length(col)==1) col<-rep(col,2)
    if (linetype=="surv" | linetype=="both"){
        segments(x$L,x$surv,x$R,x$surv,
            lty=lty[2],col=col[2],...)
        segments(x$L[-1],x$surv[-1],x$R[-n],x$surv[-n],
            lty=lty[2],col=col[2],...)
        if (mark.time & length(x$cens)>0){
            out<-StCI(x,x$cens)
            points(out$time,out$survival,pch=3,col=col[2])
        }
    }
}

summary.kmci<-function(object,...){
    ## since summary method in stats uses object, use 
    ## that, change to x because thats what I had originally 
    x<-object
    out<-data.frame(x$time,x$surv,x$lower,x$upper)
    dimnames(out)[[2]]<-c("time","survival",
           paste("lower ",100*x$conf.level,"% CL",sep=""),
           paste("upper ",100*x$conf.level,"% CL",sep=""))
    #print(out,row.names=FALSE)
    #invisible(out)
    out
}

summary.kmciLR<-function(object,...){
    x<-object
    if (is.null(x$Lin)) x$Lin<-rep(FALSE,length(x$L))
    if (is.null(x$Rin)) x$Rin<-rep(TRUE,length(x$R))
    Interval<-intChar(x$L,x$R,x$Lin,x$Rin)
    out<-data.frame(Interval,x$surv,x$lower,x$upper)
    dimnames(out)[[2]]<-c("time interval","survival",
           paste("lower ",100*x$conf.level,"% CL",sep=""),
           paste("upper ",100*x$conf.level,"% CL",sep=""))
    #print(out,row.names=FALSE)
    #invisible(out)
    out
}
#summary(norm)

summary.kmciLRgroup <- function(object, ...) {
  for (i in 1:length(object)) {
    assign(paste("ddat",i,sep="_"), cbind(summary(object[[i]]), rep(names(object)[i], length(object[[i]]$Interval))))
  }
  df <- mget(ls(pattern = "ddat"))
  out <- as.data.frame(NULL)
  for (i in 1:length(df)) {
    out <- rbind(out, df[[i]])
  }
  names(out)[5] <- object[[1]]$groupVarName
  return(out)
}

summary.kmciLRtidy <- function(object, ...) {
  if (length(object) == 1) {
    out <- summary(object[[1]])
  }
  else {
    class(object) <- "kmciLRgroup"
    out <- summary.kmciLRgroup(object, ...)
  }
  return(out)
}

StCI<-function(x,tstar,afterMax="continue",...){
   UseMethod("StCI")
}
StCI.default<-function(x,tstar,afterMax="continue",...){
    ### get survival and confidence interval at t from a survfit or kmci object
    #if (class(x)=="survfit"){
    if (is(x,"survfit")){
        x$conf.level<-x$conf.int
    }
    if (length(x$strata)>1){
     stop("does not work for objects with more than one strata")
    }
    time<-x$time
    k<-length(time)
    index<-1:k
    nt<-length(tstar)
    ### important to set I to 0, if any tstar< min(time) stays
    ### equal to 0, then we will fix them at the end
    I<-rep(0,nt)
    for (j in 1:nt){
        J<- time<=tstar[j]
        if (any(J)){
            I[j]<-max(index[J])
        }
    }
    ### afterMax determines what to do after the maximum time
    ### afterMax="continue"
    ###    - surv, lower, and upper continue at value 
    ###      at time[nt]
    ### afterMax="zero"
    ###    - surv, lower go to zero, upper continues at value 
    ###       at time[nt]
    ### afterMax="zeroNoNA"
    ###    - surv, lower go to zero, upper continues at value 
    ###       at time[nt] (unless it is NA, then take the 
    ###       last non-missing value
    ### afterMax="half"
    ###    - surv goes to half value at time[nt]
    ###    - lower goes to zero, upper continues at value
    ###       at time[nt]
    if (afterMax!="continue" && any(I==k)){
        ### default is to continue, 
        ### no need to do anything if afterMax="continue"
        if (afterMax=="zero"){
            x$surv[k]<-0
            x$lower[k]<-0
        } else if (afterMax=="zeroNoNA"){
            if (is.na(x$lower[k])) x$lower[k]<-0
            if (is.na(x$upper[k])){
                x$upper[k]<-x$upper[max(index[!is.na(x$upper)])]
            }
        } else if (afterMax=="half"){
            x$surv[k]<- .5*x$surv[k] 
            x$lower[k]<-0
        } else stop("afterMax must be 'continue',
                     'zero','zeroNoNA', or 'half' ")
    }
    ### I==0 are when tstar[j]< min(time), set all to 1
    S<-L<-U<-rep(1,nt)
    ## when I>0, i.e., tstar[j]>=min(time), plug in values
    ## note: I<-c(0,2,4), x<-1:4, then x[I] gives c(2,4), 
    ##     zeros ignored
    S[I>0]<-x$surv[I]
    L[I>0]<-x$lower[I]
    U[I>0]<-x$upper[I]

    out<-data.frame(time=tstar,survival=S,lower=L,upper=U)
    ## changing the name of the column in the data.frame 
    ## after the conf.level, was a bad idea.
    ## It is cleaner to add conf.level as an attribute. 
    #dimnames(out)[[2]]<-c("time","survival",
    #       paste("lower ",100*x$conf.level,"% CL",sep=""),
    #       paste("upper ",100*x$conf.level,"% CL",sep=""))
    attr(out,"conf.level")<- x$conf.level
    #print(out,row.names=FALSE)
    #invisible(out)
    out  
}

StCI.kmciLR<-function(x,tstar,...){
    ### get survival and confidence interval at t 
    ### from kmciLR object
    nt<-length(tstar)
    I<-rep(NA,nt)
    index<-1:length(x$surv)
    ## picki gives TRUE/FALSE vector, TRUE where tval fits 
    ## into interval
    picki<-function(tval){
        (x$L<tval & x$R>tval) | (x$L==tval & x$Lin) | 
        (x$R==tval & x$Rin)
    }
    for (j in 1:nt){
        I[j]<-index[picki(tstar[j])]
    }
    out<-data.frame(time=tstar,survival=x$surv[I],
        lower=x$lower[I],upper=x$upper[I])
    ## changing the name of the column in the data.frame 
    ## after the conf.level, was a bad idea.
    ## It is cleaner to add conf.level as an attribute. 
    #dimnames(out)[[2]]<-c("time","survival",
    #       paste("lower ",100*x$conf.level,"% CL",sep=""),
    #       paste("upper ",100*x$conf.level,"% CL",sep=""))
    #print(out,row.names=FALSE)
    attr(out,"conf.level")<- x$conf.level
    #invisible(out)
    out     
}

quantile.kmciLR<-function(x,probs=c(.25,.5,.75),...){
    lower<-upper<-surv<-rep(NA,length(probs))

    k<-length(x$surv)    
    q1<-function(x,p){
        lower.i<-min((1:k)[x$lower<=p])
        upper.i<-max((1:k)[x$upper>=p])
       if (any(x$surv==p)){
            q<-min(c(x$L[x$surv==p],x$R[x$surv==p]))
        } else if (any(x$surv<p)){
            q<- x$L[min((1:k)[x$surv<p])]            
        } else {
            q<-NA
        }
        lower<-x$L[lower.i]
        upper<-x$R[upper.i]
        c(q,lower,upper)
    }
    out<-matrix(NA,length(probs),4,dimnames=list(NULL,
        c("S(q)","q","lower","upper")))
    for (i in 1:length(probs)){
        out[i,2:4]<-q1(x,probs[i])
    }
    out[,1]<-probs
    out
}



quantile.kmci<-function(x,probs=c(.25,.5,.75),...){
    k<-length(x$surv)
    newx<-list(
        surv=c(1,x$surv),
        L=c(0,x$time),
        Lin=c(FALSE,rep(TRUE,k)),
        R=c(x$time,Inf),
        Rin=c(rep(FALSE,k+1)),
        lower=c(x$lower[1],x$lower),
        upper=c(1,x$upper))
    quantile.kmciLR(newx,probs)
}

quantile.kmciLRgroup <- function(x,probs=c(.25,.5,.75),...) {
  out <- list()
  for (i in 1:length(x)) {
    out[[i]] <- quantile(x[[i]], probs)
    names(out)[i] <- names(x)[i]
  }
  return(out)
}

quantile.kmciLRtidy <- function(x,probs=c(.25,.5,.75),...) {
  if (length(x) == 1) {
    out <- quantile(x[[1]], probs)
  }
  else {
    class(x) <- "kmciLRgroup"
    out <- quantile.kmciLRgroup(x, probs)
  }
  return(out)
}


median.kmciLR<-function(x,...){
    quantile.kmciLR(x,probs=.5)
}
median.kmci<-function(x,...){
    quantile.kmci(x,probs=.5)
}


median.kmciLRtidy <- function(x, ...) {
  quantile.kmciLRtidy(x,probs=.5)
}

median.kmciLRgroup <- function(x, ...) {
  quantile.kmciLRgroup(x,probs=.5)
}


abmm<-function(a1,b1,a2=NULL,b2=NULL){
    ## Jan 6, 2016: add NULL defaults, so can leave off a2 and b2
    ## use notation from Fan,1991
    a<-c(a1,a2)
    b<-c(b1,b2)
    ## March 25, 2016: problem if a=4 b=0, gives NaN
    ## June 14, 2016: problem if a=c(201,200) and b=c(0,0)
    ##      output 2 dimensional vector. Should be list(a=200,b=0)
    ##      also  list(a=201,b=0) would give same answers
    if (length(a)==1 | all(b==0) ){
        out<-list(a=a[length(a)],b=b[length(a)])
    }  else if (any(a==0)) {
        out<-list(a=0,b=b[1])
    }  else {
        S<-prod( a/(a+b) )
        T<-prod( (a*(a+1))/((a+b)*(a+b+1)) )
        newa<- ((S-T)*S)/(T-S^2)
        newb<- ((S-T)*(1-S))/(T-S^2)
        out<-list(a=newa,b=newb)
    }
    out
}

bpcp.mm<-function(x, alpha=0.05){
    # Jan 6, 2016: totally rewrote function according to Discrete notes. New convention!
    h<- length(x$ni)
    A<- x$ni-x$di+1
    B<- x$di

    # Calculate A+,A- and B+,B- for all unique time points on the input data set
    # (in notes: g2,g4,...,g2h)
    #  where W^+(g_{2j}) ~ Beta(A+[j],B+[j])
    #  and   W^-(g_{2j}) ~ Beta(A+[j],B+[j])
    
    Aplus<-Bplus<-Aminus<-Bminus<-rep(NA,h)

    for (i in 1:h){
        ab<-abmm(A[1:i],B[1:i])
        Aminus[i]<-ab$a
        Bminus[i]<-ab$b    
        if (i<h){
            ab<-abmm(A[1:i],B[1:i],x$ni[i+1],1)
            Aplus[i]<-ab$a
            Bplus[i]<-ab$b
        } else {
            # after last time, Wplus is a point mass at 0
            Aplus[i]<-0
            Bplus[i]<-1
        }
    }
    # There are 2h+1 intervals
    # [g0,g1), [g1,g2),...,[g_{2h}, g_{2h+1})
    #
    # In [g_j, g_{j+1}), for the upper limit we use 
    #        W^-(g_j).
    # So for the W^-(g), we need to evaluate at 
    #    g0,g1,g2,g3,g4,....,g2h
    #    But W^-(g3) = W^-(g2),
    #    because di=0 and ci=0 for (g2,g3]
    #    and similarly for all odd values
    aupper<- c(1,1,rep(Aminus[-h],each=2),Aminus[h])
    bupper<- c(0,0,rep(Bminus[-h],each=2),Bminus[h])
    upper<- qqbeta(1-alpha/2,aupper,bupper)

    # In [g_j, g_{j+1}), for the lower limit we use 
    #        W^+(g_{j+1}).
    # So for the W^+(g), we need to evaluate at 
    #    g1,g2,g3,g4,....,g_{2h+1}
    #    But W^+(g3) = W^+(g2),
    #    because di=0 and ci=0 for (g2,g3]
    #    and similarly for all odd values
    alower<- c(x$ni[1],rep(Aplus,each=2))
    blower<- c(1,rep(Bplus,each=2))
    lower<- qqbeta(alpha/2,alower,blower)

    list(upper=upper,lower=lower,alower=alower,
         blower=blower,aupper=aupper,bupper=bupper)
}




bpcpMidp.mm<-function(x,alpha=0.05,midptol=.Machine$double.eps^0.25){
    ## first calculate the usual bpcp.mm
    z<- bpcp.mm(x,alpha=alpha)
    ## extract a and b Beta parameters for lower and upper
    a1<-z$alower
    b1<-z$blower
    a2<-z$aupper
    b2<-z$bupper     

    m<-length(a1)
    lowerRootFunc<-function(x,i){
        qqbeta(alpha/2-x,a1[i],b1[i]) - 
            qqbeta(alpha/2+x,a2[i],b2[i])
    }
    
    upperRootFunc<-function(x,i){
        qqbeta(1-alpha/2-x,a1[i],b1[i]) - 
            qqbeta(1-alpha/2+x,a2[i],b2[i])
    }


     lower<-upper<-rep(NA,m)

     for (i in 1:m){
         
         if (b2[i]==0){
             # Recall:    W=U*Bl + (1-U)*Bu, where 
             #   U ~ Bernoulli(.5)
             #   Bl~ Random Variable for lower
             #   Bu~ Random variable for upper 
             # if b2[i]==0, then Bu is a point mass at 1
             # Let q(x,W) be the xth quantile of a RV W
             # so the quantile of W at alpha/2 is  q(alpha/2,W) 
             #  and 
             #  q(alpha/2,W) = q(alpha,Bl)   for all 0<alpha<1 
             lower[i]<-qqbeta(alpha,a1[i],b1[i])
             upper[i]<-1
         } else if (a1[i]==0){
             # if a1[i]==0, then Bl is a point mass at 0
             #  and 
             #  q(1-alpha/2,W) = q(1-alpha,Bu)   
             #               for all 0<alpha<1 
             lower[i]<-0
             upper[i]<-qqbeta(1-alpha,a2[i],b2[i])
         } else if (i>1 & (a1[i]==a1[i-1] & b1[i]==b1[i-1] & 
                    a2[i]==a2[i-1] & b2[i]==b2[i-1])){
             ## this if condition is just for 
             ## saving computation time
             lower[i]<-lower[i-1]
             upper[i]<-upper[i-1]
         } else {
             lowerRoot<-uniroot(lowerRootFunc,interval=
                    c(-alpha/2,alpha/2),i=i,tol=midptol)$root
             # see paper, we want 
             #    Q(alpha1,a1,b1)=Q(alpha2,ab,b2), 
             # where alpha1+alpha2=alpha
             # so we solve that using uniroot, 
             # then Pr[W<=Q]=alpha/2
             lower[i]<-qqbeta(alpha/2-lowerRoot,a1[i],b1[i])
             upperRoot<-uniroot(upperRootFunc,interval=
                   c(-alpha/2,alpha/2),i=i,tol=midptol)$root
             upper[i]<-qqbeta(1-alpha/2-upperRoot,a1[i],b1[i])
         }
    
     }

     out<-list(lower=lower,
               upper=upper,
               alower=a1,aupper=a2,blower=b1,bupper=b2)
     out  
}


#outmm<-bpcpMidp.mm(x)
#out<-bpcpMidp.mc(x,nmc=10^5)
#max( abs(outmm$lower-out$lower) )
#max( abs(outmm$upper-out$upper) )


bpcpControl<-function(midpMMTol=.Machine$double.eps^0.25, 
    seed=49911,tolerance=.Machine$double.eps^0.5){
    # if you put seed=NA change it to seed=NULL
    if (!is.null(seed) & is.na(seed)){ seed<-NULL } 
    if (tolerance<=0){ 
        stop("tolerance must be positive. 
             It is the lowest positive value such that if 
             abs(x-y) is less than tolerance, 
             then numerics x and y are treated as equal")  }
    list(midpMMTol=midpMMTol,seed=seed,tolerance=tolerance)
}


# bpcp.mc is a function to calculate the bpcp CIs by Monte Carlo simulation
bpcp.mc<-function(x,nmc=100,alpha=.05, testtime=0, DELTA=0, midp=FALSE){
    # Jan 6, 2016: totally rewrote function according to Discrete notes. New convention!
    ### for each time t_j there are 2 intervals
    ### representing
    ## [t_{j-1},t_j-Delta)
    ## [t_j-Delta, t_j)
    ##
    ## and at the end add [t_k, Inf)
    ## 
    ## if Delta=0 then the second interval of the pair is 
    ## not needed
    ## but we keep it in and delete in bpcp after the 
    ## call to this .calc function 
    k<-length(x$time)
    lower<-upper<-rep(1,2*k+1)
    S<-rep(1,nmc)
    q<-function(slo,shi,Midp=midp){
        if (Midp){
            S<-c(slo,shi)
            quantile(S,probs=c(alpha/2,1-alpha/2))
        } else {
            c( quantile(slo, probs=alpha/2), quantile(shi,1-alpha/2) )
        }
    }

    # function to pick out time points
    # so t_(j-1) = t_{j-1}, for j=1,..k
    t_<-function(j){
        tt<-c(0,x$time)
        tt[j+1]
    }
    Smc<-list(Slo=NA,Shi=NA)
    for (j in 1:k){
        ## case 1: t_j is a death time (and perhaps censor 
        ## time also)
        if (x$di[j]>0){
            Shi<-S
            ## Slo=look ahead one failure
            Slo<-S*rbeta(nmc,x$ni[j],1)
            lohi<-q(Slo,Shi)
            ## for  (t_{j-1},t_j-Delta]
            if (t_(j-1)< testtime & 
               testtime<=t_(j)-DELTA) Smc<-list(Slo=Slo,Shi=Shi)
            lower[2*j-1]<-lohi[1]
            upper[2*j-1]<-lohi[2]
            Slo<-S*rbeta(nmc,x$ni[j] - x$di[j] + 1,x$di[j])
            lohi<-q(Slo,Shi)
            ## for (t_j-Delta, t_j]
            if (t_(j)-DELTA< testtime & 
               testtime<=t_(j)) Smc<-list(Slo=Slo,Shi=Shi)
            lower[2*j]<-lohi[1]
            upper[2*j]<-lohi[2]
            S<-Slo
            Shi<-Slo
            lohi<-q(Slo,Shi)
         }  else {
            Shi<-S
            ## Slo=look ahead one failure
            Slo<-S*rbeta(nmc,x$ni[j],1)
            lohi<-q(Slo,Shi)
            ## for (t_{j-1},t_j]
            if (t_(j-1)< testtime & 
                  testtime<=t_(j)) Smc<-list(Slo=Slo,Shi=Shi)
            lower[(2*j-1):(2*j)]<-lohi[1]
            upper[(2*j-1):(2*j)]<-lohi[2]
        }  
    }
    ## for (t_k, Inf)
    if (testtime>t_(k)) Smc<-list(Slo=S,Shi=rep(0,nmc))
    lohi<- q(rep(0,nmc),S)
    upper[2*k+1]<-lohi[2]
    lower[2*k+1]<-lohi[1]
    list(upper=upper,lower=lower, Smc=Smc)
}






bpcp<-function(time,status,nmc=0,alpha=.05,Delta=0,stype="km", midp=FALSE, monotonic=NULL, control=bpcpControl()){
    # Jan 6, 2016: totally rewrote function according to Discrete notes. New convention!
    if (is.null(monotonic)){
        if (nmc==0){
            monotonic<-TRUE
        } else {
            monotonic<-FALSE
        }
    }
    midpMMTol<-control$midpMMTol
    # so we get the same answer with the same data set, 
    # default to set seed...control$seed=NULL is for
    # simulations on Monte Carlo
    if (nmc>0 & !is.null(control$seed)) set.seed(control$seed)
    ##
    ## get Kaplan-Meier and Greenwood variances
    x<-kmgw.calc(time,status,keepCens=TRUE)
    k<-length(x$time)
    # Check Delta:
    if (Delta<0) stop("Delta must be greater than 0")
    minTimeDiff<-min(diff(c(0,x$time)))
    # instead of using minTimeDiff<Delta, use the following
    # in the if condition, so that machine error does not 
    # create problems (see default for all.equal tolerance)
    tolerance<-control$tolerance
    if (minTimeDiff-Delta+
         tolerance<0) stop(
              "Either negative times or 
              Delta is not less than or equal to 
              the minimum difference in times")
    ## each time, t_j, represents 2 intervals
    ## [t_{j-1},t_j-Delta)
    ## [t_j-Delta, t_j)
    ##  add on KM at t=0,
    ##  for the two intervals: [0,t1-Delta), (t1-Delta,t1]
    KM<-c(rep(1,2),rep(x$KM[-k],each=2),x$KM[k])
    L<-c(0,rep(x$time,each=2)) - c(0,rep(c(Delta,0),k))
    R<-c(rep(x$time,each=2),Inf) - c(rep(c(Delta,0),k),0)
    Lin<- rep(TRUE,2*k+1)
    Rin<- rep(FALSE,2*k+1)

    if (midp){
        if (nmc==0){ 
           hilo<-bpcpMidp.mm(x,alpha=alpha,midptol=midpMMTol)
        } else hilo<-bpcp.mc(x,nmc,alpha,midp=TRUE)
    } else {
        if (nmc==0){ hilo<-bpcp.mm(x,alpha=alpha)
        } else hilo<-bpcp.mc(x,nmc,alpha)
     }

    ## we do not need to keep all the intervals in all cases
    ## 1) if Delta==0 then do not need 2nd interval 
    ## of each pair
    if (Delta==0){ 
        keep<-c(rep(c(TRUE,FALSE),k),TRUE)
    } else {
    ## 2) if Delta>0, sometimes you have deaths or 
    ##    censoring in adjascent 
    ##    intervals. For example if Delta=1 and you have 
    ##    deaths in 
    ##    (0,1]  then the death "time" is 1 
    ##    and the 2 intervals for t_j=1 are
    ##        [t0,t1-Delta) = [0,0)
    ##        [t1-Delta, t1)    =[0,1)
    ##    and the first interval is not needed. If the 
    ##    next death or censoring 
    ##    time was t2=3
    ##    the next 2 intervals are:
    ##                  [t1,t2-Delta) = [1,3-1) =[1,2)
    ##                  [t2-Delta, t2) = [2,3)
    ##    and we do not need to delete any. 
    ##    
    ##     So basically if t_j-Delta=t_{j-1} then we 
    ##     delete the first interval 
    ##     of the pair.
    #    We cannot use the following code for keepfirst
    #    because there may be machine error problems
    #    such as 3.1-.1 != 3   being TRUE 
    #        keepfirst<- diff(c(0,x$time))!=Delta
    #   So we call numerics within .Machine$double.eps^0.5 
    #   as equal   (see default for all.equal)
    # 
    keepfirst<- abs(diff(c(0,x$time))-
                rep(Delta,length(x$time)))>=
                            tolerance
    keep<-rep(TRUE,2*k+1)
    first<-c(rep(c(TRUE,FALSE),k),FALSE)

    keep[first]<-keepfirst
    }
    if (stype=="km"){ 
        SURV<-KM[keep]
    } else {    
        if (stype!="mue"){
            warning("assuming stype='mue' ")
            esimate<-"mue"
        }
        ## use recursive call... fix this later to 
        ## make it faster
        outtemp<-bpcp(time,status,nmc,
            alpha=1,Delta,stype="km",midp=midp)
        SURV<- .5*outtemp$lower + .5*outtemp$upper
    }
    ## create list of beta parameters to go with the 
    ## lower and upper CIs
    betaParms<-NULL
    if (nmc==0) betaParms<-list(
        alower=hilo$alower[keep],
        blower=hilo$blower[keep],
        aupper=hilo$aupper[keep],
        bupper=hilo$bupper[keep])
    lower<-hilo$lower[keep]
    upper<-hilo$upper[keep]
    if (monotonic){
        lower<-cummin(lower)
        upper<-cummin(upper)
    }
    out<-list(cens=getmarks(time,status),surv=SURV,
        lower=lower,
        upper=upper,
        L=L[keep],Lin=Lin[keep],R=R[keep],Rin=Rin[keep],
        Interval=intChar(L,R,Lin,Rin)[keep],stype=stype, 
        betaParms=betaParms, conf.level=1-alpha)
    class(out)<-"kmciLR"
    out
}

#Function to create "tidy" output from bpcp function
#Transform kmciLR object into dataframe 
#Takes in a list of kmciLR objects
tidykmciLR <- function(x) {
  #Create empty dataframe to store output
  tidyout <- data.frame(NULL)
  
  #For each group, "tidy" the output into a new dataframe
  if (attr(x, "class") %in% c("kmciLRtidy", "kmciLRgroup")) {
    num <- length(x)
    for (i in 1:length(x)) {
      new <- with(x[[i]], data.frame(time = sort(c(L, R)), surv = rep(surv, each = 2), 
                                     lower = rep(lower, each = 2), upper = rep(upper, each = 2)))
      #Add group variable if it exists
      if (num != 1) {
        new$group <- names(x[i])
      }
      #Add to tidy dataframe
      tidyout <- rbind(tidyout, new)
    }
  }
  else if (attr(x, "class") == "kmciLR") {
    tidyout <- with(x, data.frame(time = sort(c(L, R)), surv = rep(surv, each = 2), 
                                  lower = rep(lower, each = 2), upper = rep(upper, each = 2)))
  }
  return(tidyout)
}

bpcpfit <- function(time, ...) {
  UseMethod("bpcpfit")
}

#Takes same inputs as bpcp function
bpcpfit.formula <- function(formula, data, subset, na.action, ...) {
  #Borrowed from survfit
  Call <- match.call()
  Call[[1]] <- as.name('bpcpfit')  
  indx <- match(c('formula', 'data'), names(Call), nomatch=0)
  #Make sure formula is given
  if (indx[1]==0) {
    stop("a formula argument is required")
  }
  if (missing(data)) 
    data <- environment(formula)
  #Change to model.frame format using data
  temp <- model.frame(formula, data=data)
  #Evaluate
  mf <- eval.parent(temp)
  #Get terms of formula
  Terms <- terms(formula)
  #Make sure there are no interaction terms present in the formula
  ord <- attr(Terms, 'order')
  if (length(ord) & any(ord !=1)) {
    stop("Interaction terms are not valid for this function")
  }
  n <- nrow(mf)
  #Get out Y variable (time and censoring)
  Y <- model.extract(mf, 'response')
  #Ensure a Surv object is provided to the formula
  if (!is.Surv(Y)) {
    stop("Response must be a survival (Surv) object")
  }
  #Get labels from formula (group variable)
  ll <- attr(Terms, 'term.labels')
  if (length(ll) > 1) {
    stop("Only a treatment/grouping variable can be specified in this function. No other covariates should be included.")
  }
  #Determine if group variable was provided
  if (length(ll) == 0) {
    output <- do.call("bpcpfit.default", list(Y[ ,1], Y[ ,2], ...))
  }
  else {
    X <- strata(mf[ll], shortlabel = TRUE)
    #Combine outcome and response variables
    Z <- data.frame(cbind(Y, X))
    #Split my levels of group (treatment) variable
    newZ <- split(Z, Z$X)
    
    #Iterate through each level of the treatment variable, and do the bpcp function of the corresponding data. Store the results in a list
    for (i in 1:length(newZ)) {
      if (length(unique(newZ[[i]]$status)) > 2) {
        stop("Interval censoring is not supported by bpcp or bpcpfit.")
      }
    }
    output <- do.call("bpcpfit.default", list(Z[ ,1], Z[ ,2], Z[ ,3], ...))
    names(output) <- levels(factor(X))
    for (i in 1:length(output)) {
      output[[i]]$groupVarName <- ll
    }
  }
  return(output)
}

bpcpfit.default <- function(time, status = NULL, group = NULL, formula=NULL, nmc=0, alpha=NULL, conf.level=0.95, Delta=0, stype="km", midp=FALSE,
                            monotonic=NULL, control=bpcpControl(), plotstyle = "ggplot", data=NULL, subset=NULL, na.action=NULL, ...) {
  if (is.null(time)) {
    stop("Time is required.")
  }
  if (!is.null(alpha)) {
    print("Warning: alpha is out of date. Use conf.level instead. Setting conf.level to 1-alpha and ignoring conf.level.")
  }
  else {
    alpha <- 1-conf.level
  }
  Call <- match.call()
  if (is.null(group)) {
    if (length(unique(status)) > 2) {
      stop("Interval censoring is not supported by bpcp or bpcpfit.")
    }
    out <- bpcp(time, status, nmc=nmc, alpha=alpha, Delta=Delta, stype=stype, midp=midp,
                monotonic=monotonic, control=control)
    out$num <- length(time)
    out$events <- length(which(status == max(status)))
    if (plotstyle == "ggplot") {
      results <- list()
      results <- append(results, list(out))
      class(results) <- "kmciLRtidy"
    }
    else {
      results <- out
    }
  }
  else {
    results <- list()
    Z <- as.data.frame(cbind(time, status, group))
    names(Z) <- c("V1", "V2", "V3")
    newZ <- split(Z, Z$V3)
    results <- list()
    #Iterate through each level of the treatment variable, and do the bpcp function of the corresponding data. Store the results in a list
    for (i in 1:length(newZ)) {
      if (length(unique(newZ[[i]]$V2)) > 2) {
        stop("Interval censoring is not supported by bpcp or bpcpfit.")
      }
      results <- append(results, list(bpcp(newZ[[i]]$V1, newZ[[i]]$V2, nmc=nmc, alpha=alpha, Delta=Delta, stype=stype, midp=midp,
                                           monotonic=monotonic, control=control)))
      results[[i]]$num <- length(newZ[[i]]$V1)
      results[[i]]$events <- length(which(newZ[[i]]$V2 == max(newZ[[i]]$V2)))
      results[[i]]$groupVarName <- rev(strsplit(as.character(Call["group"]), "$", fixed = TRUE)[[1]])[1]
    }
    #Name each output from the bpcp function with the corresponding treatment
    names(results) <- levels(factor(group))
    if (plotstyle == "ggplot") {
      class(results) <- "kmciLRtidy"
    }
    else if (plotstyle == "standard") {
      class(results) <- "kmciLRgroup"
    }
    else {
      stop('Plot style must be either "ggplot" or "standard"')
    }
  }
  return(results)
}


print.kmciLRgroup <- function(x, ...) {
  #Store results in a matrix
  out<-matrix(NA,length(x),5,dimnames=list(NULL, c("n", "events", "median", paste0(x[[1]]$conf.level, "LCL"), paste0(x[[1]]$conf.level, "UCL"))))
  #For each level of treatment variable, also get number of subjects and number of events and print those as well
  for (i in 1:length(x)) {
    out[i,3:5]<- median(x)[[i]][2:4]
    out[i, 1] <- x[[i]]$num
    out[i, 2] <- x[[i]]$events
  }
  #Row names are levels of treatment variable
  rownames(out) <- names(x)
  invisible(x)
  print(out)
}

#Output mirrors output of survfit function
print.kmciLRtidy <- function(x, ...) {
  if (length(x) == 1) {
    print(x[[1]])
  }
  else{
    class(x) <- "kmciLRgroup"
    print.kmciLRgroup(x)
  }
  invisible(x)
}

print.kmciLR <- function(x, ...) {
  #Store results in a matrix
  out<-matrix(NA,1,5,dimnames=list(NULL, c("n", "events", "median", paste0(x$conf.level, "LCL"), paste0(x$conf.level, "UCL"))))
  out[1,3:5]<- median(x)[2:4]
  out[1, 1] <- x$num
  out[1, 2] <- x$events
  
  #Row names are levels of treatment variable
  row.names(out) <- NULL
  invisible(x)
  print(out)
}


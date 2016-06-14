## Thanks to Marcel Wolbers for finding this problem.
## It used to produce a warning, but the answers were correct.
## Now should not produce a warning

library(bpcp)
library(survival)



sim.exp <- function(n,prob1,prob2,cens.rate){
  # prob1 and prob2 are the absolute risks for an event at time t=1
  rate1 <- -log(1-prob1); rate2 <- -log(1-prob2)
  n1 <- n%/%2+n%%2; n2 <- n%/%2
  ttev.uncens <- rexp(n,rate=rep(c(rate1,rate2),times=c(n1,n2)))
  ttcens <- rexp(n,rate=cens.rate)  
  data.frame(arm=rep(c(0,1),c(n1,n2)),ttev=pmin(ttev.uncens,ttcens),ev=ifelse(ttev.uncens <=ttcens,1,0))
}

set.seed(1)
d <- sim.exp(400,prob1=0.1,prob2=0.05,cens.rate=0.05)
summary(survfit(Surv(ttev,ev)~arm,data=d),time=1)
bpcp2samp(d$ttev,d$ev,d$arm,testtime=1,parmtype="difference")

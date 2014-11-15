#######################################################
# R script for 11.5 A Robust Regression Model
#######################################################

require(arm)
require(LearnBayes)

data(election)
attach(election)
y=sqrt(buchanan)
x=sqrt(perot)
N=length(y)

data=list("N","y","x")
inits = function() {list(b=c(0,0),tau=1)}
parameters <- c("tau","lam","b")
robust.sim <- bugs (data, inits, parameters, "robust.bug", n.chains=3, n.iter=1000)
print(robust.sim)

attach.bugs(robust.sim)
xo=seq(18,196,2)
X0=cbind(1,xo)
meanresponse=b%*%t(X0)
meanp=apply(meanresponse,2,quantile,c(.05,.5,.95))
plot(sqrt(perot),sqrt(buchanan))
lines(xo,meanp[2,])
lines(xo,meanp[1,],lty=2)
lines(xo,meanp[3,],lty=2)

#######################################################


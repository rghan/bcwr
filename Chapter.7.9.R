#################################################
# Section 7.9  Bayesian Sensitivity Analysis
#################################################

 library(LearnBayes)
 data(hearttransplants)
 attach(hearttransplants)

 datapar = list(data = hearttransplants, z0 = 0.53)

 start = c(4, -7)
 fitgibbs = gibbs(poissgamexch, start, 1000, c(1,.15), datapar)


sir.old.new=function(theta, prior, prior.new)
{
log.g=log(prior(theta))
log.g.new=log(prior.new(theta))
wt=exp(log.g.new-log.g-max(log.g.new-log.g))
probs=wt/sum(wt)
n=length(probs)
indices=sample(1:n,size=n,prob=probs,replace=TRUE)
theta[indices]
}

prior=function(theta)
0.53*exp(theta)/(exp(theta)+0.53)^2
prior.new=function(theta)
5*exp(theta)/(exp(theta)+5)^2

log.alpha=fitgibbs$par[, 1]
log.alpha.new=sir.old.new(log.alpha, prior, prior.new)

############ drawing figure

library(lattice)
draw.graph=function()
{
LOG.ALPHA=data.frame("prior",log.alpha)
names(LOG.ALPHA)=c("Prior","log.alpha")
LOG.ALPHA.NEW=data.frame("new.prior",log.alpha.new)
names(LOG.ALPHA.NEW)=c("Prior","log.alpha")
D=densityplot(~log.alpha,group=Prior,data=rbind(LOG.ALPHA,LOG.ALPHA.NEW),
   plot.points=FALSE,main="Original Prior and Posterior (solid), 
   New Prior and Posterior (dashed)",
    lwd=4,adjust=2,lty=c(1,2),xlab="log alpha",xlim=c(-3,5),col="black")
update(D, panel=function(...){
    panel.curve(prior(x),lty=1,lwd=2,col="black")
    panel.curve(prior.new(x),lty=2, lwd=2,col="black")
    panel.densityplot(...)
})}

draw.graph()

############# does posterior on true rates depend on prior?

prob.int.rate=function(i,log.alpha)
{
     lami = rgamma(1000, y[i] + exp(log.alpha), e[i] + exp(log.alpha)/mu)
     quantile(lami, c(0.05, 0.95))
}

ind=c(1,25,57,89)

ci=sapply(ind,prob.int.rate,log.alpha)
matplot(matrix(log(e[ind]),2,length(ind),byrow=T),ci,type="l",
   lwd=4,lty=1,xlim=c(6,9.7),ylim=c(0,.002),ylab="True Rate",xlab="log(e)")

ci=sapply(ind,prob.int.rate,log.alpha.new)
matplot(matrix(log(e[ind]),2,length(ind),byrow=T)+.05,ci,type="l",lwd=4,lty=3,
  add=TRUE)


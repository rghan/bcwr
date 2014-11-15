#######################################################################
#  R script for Section 11.6  Estimating Career Trajectories
#######################################################################

require(arm)
require(LearnBayes)

data(sluggerdata)
s=careertraj.setup(sluggerdata)
N=s$N; S=s$T; y=s$y; n=s$n; x=s$x

mean = c(0, 0, 0)
Omega=diag(c(.1,.1,.1))
prec=diag(c(1.0E-6,1.0E-6,1.0E-6))

beta0=matrix(c(-7.69,.350,-.0058),nrow=10,ncol=3,byrow=TRUE)
mu.beta0=c(-7.69,.350,-.0058)
R0=diag(c(.1,.1,.1))

data=list("N","S","y","n","x","mean","Omega","prec")
inits = function() {list(beta=beta0,mu.beta=mu.beta0,R=R0)}
parameters <- c("beta")
career.sim <- bugs (data, inits, parameters, "career.bug", 
   n.chains=1, n.iter=10000, n.thin=1, codaPkg=TRUE)

career.coda = read.bugs(career.sim)
windows(record=TRUE)
plot(career.coda,ask=TRUE)
summary(career.coda)
densityplot(career.coda)
###########################################################################

career.sim <- bugs (data, inits, parameters, "career.bug", 
   n.chains=1, n.iter=50000, n.thin=1)

peak.age=matrix(0,50000,10)
for (i in 1:10)
peak.age[,i]=-career.sim$sims.list$beta[,i,2]/2/
 career.sim$sims.list$beta[,i,3]

dimnames(peak.age)[[2]]=c("Aaron","Greenberg", "Killebrew", "Mantle","Mays",
      "McCovey" ,"Ott", "Ruth", "Schmidt", "Sosa") 
densityplot(as.mcmc(peak.age),plot.points=FALSE,col="black",lwd=2)

summary(as.mcmc(peak.age))

Quantiles for each variable:

          2.5%  25%  50%  75%  98%
Aaron     31.3 32.2 32.8 33.4 34.8
Greenberg 29.7 31.2 32.1 33.2 35.8
Killebrew 26.1 27.1 27.6 28.0 28.9
Mantle    27.1 27.9 28.3 28.7 29.8
Mays      28.8 29.8 30.2 30.8 31.8
McCovey   26.8 28.6 29.2 29.7 30.7
Ott       26.3 27.3 27.8 28.4 29.8
Ruth      30.1 31.0 31.5 32.0 33.1
Schmidt   27.6 28.7 29.2 29.7 30.7
Sosa      30.8 32.1 32.9 33.9 35.9




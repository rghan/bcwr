#################################################
# R script for Section 11.4 A Change-Point Model
#################################################

require(arm)

N=112
D=c(4,5,4,1,0,4,3,4,0,6,
3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,3,4,2,5,2,2,3,4,2,1,3,2,
1,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,2,
0,1,1,1,0,1,0,1,0,0,0,2,1,0,0,0,1,1,0,2,
2,3,1,1,2,1,1,1,1,2,4,2,0,0,0,1,4,0,0,0,
1,0,0,0,0,0,1,0,0,1,0,0)
data=list("N","D")
inits = function() {list(b=c(0,0),changeyear=50)}
parameters <- c("changeyear","b")
coalmining.sim <- bugs (data, inits, parameters, "coalmining.bug", 
   n.chains=3, n.iter=1000, codaPkg=TRUE)
coalmining.coda = read.bugs(coalmining.sim)
summary(coalmining.coda)
xyplot(coalmining.coda)
acfplot(coalmining.coda)
densityplot(coalmining.coda,col="black")



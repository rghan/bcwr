# Bayesian Computation with R (Second Edition)
# by Jim Albert

# 3 Single-Parameter Models

# 3.1 Introduction

# 3.2 Normal Distribution with Known Mean but Unknown Variance

data(footballscores)
attach(footballscores)
d = favorite - underdog - spread
n = length(d)
v = sum(d^2)

P = rchisq(1000, n)/v
s = sqrt(1/P)
hist(s, main ="")

# Fig. 3.1. Histogram of simulated sample of the standard deviation σ of differences between
# game outcomes and point spreads.

quantile(s, probs = c(0.025, 0.5, 0.975))

# 3.3 Estimating a Heart Transplant Mortality Rate

alpha = 16
beta = 15174
yobs = 1
ex = 66
y = 0:10
lam = alpha/beta
py = dpois(y, lam*ex)*dgamma(lam, shape = alpha, rate = beta)/
  dgamma(lam, shape= alpha + y, rate = beta + ex)
cbind(y, round(py, 3))

lambdaA = rgamma(1000, shape = alpha + yobs, rate = beta + ex)

ex = 1767
yobs = 4
y = 0:10
py = dpois(y, lam * ex) * dgamma(lam, shape = alpha, rate = beta)/
  dgamma(lam, shape = alpha + y, rate = beta + ex)
cbind(y, round(py, 3))

lambdaB = rgamma(1000, shape = alpha + yobs, rate = beta + ex)

par(mfrow = c(2, 1))
plot(density(lambdaA), main = "HOSPITAL A", xlab = "lambdaA", lwd = 3)
curve(dgamma(x, shape = alpha, rate = beta), add = TRUE)
legend("topright", legend = c("prior", "posterior"), lwd = c(1,3))
plot(density(lambdaB), main = "HOSPITAL B", xlab = "lambdaB", lwd = 3)
curve(dgamma(x, shape = alpha, rate = beta), add = TRUE)
legend("topright", legend = c("prior", "posterior"), lwd = c(1,3))

Fig. 3.2. Prior and posterior densities for heart transplant death rate for two hospitals.


# 3.4 An Illustration of Bayesian Robustness

quantile1 = list(p = 0.50, x = 100)
quantile2 = list(p = 0.95, x = 120)
normal.select(quantile1, quantile2)

mu = 100
tau = 12.16
sigma = 15
n = 4
se = sigma/sqrt(4)
ybar = c(110, 125, 140)
tau1 = 1/sqrt(1/se^2 + 1/tau^2)
mu1 = (ybar/se^2 + mu/tau^2) * tau1^2
summ1 = cbind(ybar, mu1, tau1)
summ1

tscale = 20/qt(0.95, 2)
tscale

par(mfrow = c(1,1))
curve(1/tscale*dt((x - mu)/tscale, 2), from = 60, to = 140, 
      xlab = "theta", ylab = "Prior Density") 
curve(dnorm(x, mean = mu, sd = tau), add = TRUE, lwd = 3) 
legend("topright", legend = c("t density", "normal density"), lwd = c(1, 3))

Fig. 3.3. Normal and t priors for representing prior opinion about a person’s true IQ score.

norm.t.compute = function(ybar) {
  theta = seq(60, 180, length = 500)
  like = dnorm(theta, mean = ybar, sd = sigma/sqrt(n))
  prior = dt((theta - mu)/tscale, 2)
  post = prior * like
  post = post/sum(post)
  m = sum(theta * post)
  s = sqrt(sum(theta^2 * post) - m^2)
  c(ybar, m, s) 
}
summ2 = t(sapply(c(110, 125, 140), norm.t.compute))
dimnames(summ2)[[2]] = c("ybar", "mu1 t", "tau1 t") 
summ2

cbind(summ1,summ2)

theta=seq(60, 180, length=500)
normpost = dnorm(theta, mu1[3], tau1)
normpost = normpost/sum(normpost)
plot(theta, normpost, type = "l", lwd = 3, ylab = "Posterior Density")
like = dnorm(theta, mean = 140, sd = sigma/sqrt(n))
prior = dt((theta - mu)/tscale, 2)
tpost = prior * like / sum(prior * like)
lines(theta, tpost)
legend("topright", legend = c("t prior", "normal prior"), lwd = c(1,3))

Fig. 3.4. Posterior densities for a person’s true IQ using normal and t priors for an extreme observation.


# 3.5 Mixtures of Conjugate Priors

Fig. 3.5. Mixture of beta densities prior distribution that reflects belief that a coin is biased.

probs = c(0.5, 0.5)
beta.par1 = c(6, 14)
beta.par2 = c(14, 6)
betapar = rbind(beta.par1, beta.par2)
data = c(7, 3)
post = binomial.beta.mix(probs, betapar, data)
post

curve(post$probs[1]*dbeta(x, 13, 17) + post$probs[2]*dbeta(x,21,9),
      from = 0, to = 1, lwd = 3, xlab = "P", ylab = "DENSITY")
curve(0.5*dbeta(x, 6, 12)+ 0.5*dbeta(x, 12, 6), 0, 1,add = TRUE)
legend("topleft", legend = c("Prior", "Posterior"), lwd = c(1, 3))

Fig. 3.6. Prior and posterior densities of a proportion for the biased coin example. 


# 3.6 A Bayesian Test of the Fairness of a Coin

> pbinom(5, 20, 0.5)
[1] 0.02069473

> n = 20
>y=5
> a = 10
> p = 0.5
> m1 = dbinom(y, n, p) * dbeta(p, a, a)/dbeta(p, a + y, a + n - + y)
> lambda = dbinom(y, n, p)/(dbinom(y, n, p) + m1) > lambda
[1] 0.2802215

> pbetat(p,.5,c(a,a),c(y,n-y))
$bf
[1] 0.3893163
$post
[1] 0.2802215

> prob.fair=function(log.a)
  +{
    + a = exp(log.a)
    + m2 = dbinom(y, n, p) * dbeta(p, a, a)/ + dbeta(p, a + y, a + n - y) + dbinom(y, n, p)/(dbinom(y, n, p) + m2) +}

> n = 20; y = 5; p = 0.5
> curve(prob.fair(x), from=-4, to=5, xlab="log a", > ylab="Prob(coin is fair)", lwd=2)

Fig. 3.7. Posterior probability that a coin is fair graphed against values of the prior parameter log a.

> n=20
> y=5
> a=10
> p=.5
> m2=0
> for (k in 0:y)
  +  m2=m2+dbinom(k,n,p)*dbeta(p,a,a)/dbeta(p,a+k,a+n-k)
> lambda=pbinom(y,n,p)/(pbinom(y,n,p)+m2)
> lambda
[1] 0.2184649


# 3.7 Further Reading

# 3.8 Summary of R Functions

# 3.9 Exercises

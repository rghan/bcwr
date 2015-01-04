# Bayesian Computation with R (Second Edition)
# by Jim Albert

library(LearnBayes)
library(coda)

# 6 Markov Chain Monte Carlo Methods

# 6.1 Introduction

# 6.2 Introduction to Discrete Markov Chains

P <- matrix(c(0.50, 0.50, 0.00, 0.00, 0.00, 0.00, 
              0.25, 0.50, 0.25, 0.00, 0.00, 0.00, 
              0.00, 0.25, 0.50, 0.25, 0.00, 0.00,
              0.00, 0.00, 0.25, 0.50, 0.25, 0.00,
              0.00, 0.00, 0.00, 0.25, 0.50, 0.25,
              0.00, 0.00, 0.00, 0.00, 0.50, 0.50),
            nrow = 6, ncol = 6, byrow = TRUE)
P

s <- array(0, c(50000, 1))

s[1] <- 3
for (j in 2:50000)
  s[j] <- sample(1:6, size = 1, prob = P[s[j - 1], ])

m <- c(500, 2000, 8000, 50000)
for (i in 1:4)
  print(table(s[1:m[i]])/m[i])

w <- matrix(c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1), nrow = 1, ncol = 6)
w%*%P

# 6.3 Metropolis-Hastings Algorithms

# 6.4 Gibbs Sampling

# 6.5 MCMC Output Analysis

# 6.6 A Strategy in Bayesian Computing

# 6.7 Learning About a Normal Population from Grouped Data

# Table 6.1. Grouped frequency data for heights of male students at a college.

d <- list(int.lo = c(-Inf, seq(66, 74, by = 2)),
          int.hi = c(seq(66, 74, by = 2), Inf),
          f = c(14, 30, 49, 70, 33, 15))

y <- c(rep(65, 14), rep(67, 30), rep(69, 49), rep(71, 70), rep(73, 33), rep(75, 15))
mean(y)

log(sd(y))

start <- c(70, 1)
fit <- laplace(groupeddatapost, start, d)
fit

modal.sds <- sqrt(diag(fit$var))

proposal <- list(var = fit$var, scale = 2)
fit2 <- rwmetrop(groupeddatapost, proposal, start, 10000, d)

fit2$accept

post.means <- apply(fit2$par, 2, mean)
post.sds <- apply(fit2$par, 2, sd)

cbind(c(fit$mode), modal.sds)

cbind(post.means, post.sds)

mycontour(groupeddatapost, c(69, 71, 0.6, 1.3), d, xlab = "mu", ylab = "log sigma")
points(fit2$par[5001:10000, 1], fit2$par[5001:10000, 2])

# Fig. 6.1. Contour plot of posterior of μ and log(σ) for grouped data example. A
# simulated sample of 5000 draws of the posterior is also shown.

# 6.8 Example of Output Analysis

start <- c(65,1)
proposal <- list(var = fit$var, scale = 0.2)

bayesfit <- rwmetrop(groupeddatapost, proposal, start, 10000, d)

dimnames(bayesfit$par)[[2]] <- c("mu", "log sigma")
xyplot(mcmc(bayesfit$par[-c(1:2000), ]), col = "black")

# Fig. 6.2. Trace plots of simulated draws of μ and log σ for an MCMC chain with
# poor choices for the starting value and scale factor. The first 2000 draws have been
# discarded to remove the burn-in factor.

par(mfrow = c(2, 1))
autocorr.plot(mcmc(bayesfit$par[-c(1:2000), ]), auto.layout = FALSE)

# Fig. 6.3. Autocorrelation plots of simulated draws of μ and log(σ) for an MCMC
# chain with poor choices for the starting value and scale factor.

summary(mcmc(bayesfit$par[-c(1:2000), ]))
batchSE(mcmc(bayesfit$par[-c(1:2000), ]), batchSize = 50)

start <- c(70, 1)
proposal <- list(var = fit$var, scale = 2.0)
bayesfit <- rwmetrop(groupeddatapost, proposal, start, 10000, d)

dimnames(bayesfit$par)[[2]] <- c("mu","log sigma")
sim.parameters <- mcmc(bayesfit$par[-c(1:2000), ])
xyplot(mcmc(bayesfit$par[-c(1:2000), ]), col = "black")

# Fig. 6.4. Trace plot of simulated draws of μ for an MCMC chain with good choices
# for the starting value and scale factor.

par(mfrow = c(2, 1))
autocorr.plot(sim.parameters, auto.layout = FALSE)

# Fig. 6.5. Autocorrelation plot of simulated draws of μ for an MCMC chain with
# good choices for the starting value and scale factor.

summary(sim.parameters)
batchSE(sim.parameters, batchSize = 50)

# 6.9 Modeling Data with Cauchy Errors

str(darwin)
mean(darwin$difference)
log(sd(darwin$difference))

laplace(cauchyerrorpost, c(21.6, 3.6), darwin$difference)

laplace(cauchyerrorpost, 0.1*c(21.6, 3.6), darwin$difference)$mode

c(24.70 - 4*sqrt(34.960), 24.70 + 4*sqrt(34.960))
c( 2.77 - 4*sqrt( 0.138),  2.77 + 4*sqrt( 0.138))

par(mfrow = c(1, 1))
mycontour(cauchyerrorpost, c(-10, 60, 1, 4.5), darwin$difference, xlab = "mu", ylab = "log sigma")

# Fig. 6.6. Contour plot of the posterior of μ and log σ for the Cauchy error model problem.

fitlaplace <- laplace(cauchyerrorpost, c(21.6, 3.6), darwin$difference)
mycontour(lbinorm, c(-10, 60, 1, 4.5), list(m = fitlaplace$mode, v = fitlaplace$var), xlab = "mu", ylab = "log sigma")

# Fig. 6.7. Contour plot of the normal approximation to the posterior of μ and log(σ) for the Cauchy error model problem.

proposal <- list(var = fitlaplace$var, scale = 2.5)
start <- c(20, 3)
m <- 1000
s <- rwmetrop(cauchyerrorpost, proposal, start, m, darwin$difference)
mycontour(cauchyerrorpost, c(-10, 60, 1, 4.5), darwin$difference, xlab = "mu", ylab = "log sigma")
points(s$par[, 1], s$par[, 2])

# Fig. 6.8. Contour plot of the posterior of μ and log σ with a simulated sample for the Cauchy error model problem.

# Fig. 6.9. Posterior density of μ using the normal approximation and simulated draws from the Metropolis random walk chain.

# Fig. 6.10. Posterior density of log σ using the normal approximation and simulated draws from the Metropolis random walk chain.

fitgrid <- simcontour(cauchyerrorpost, c(-10, 60, 1, 4.5), darwin$difference, 50000)
proposal <- list(var = fitlaplace$var, scale = 2.5)
start <- c(20, 3)
fitrw <- rwmetrop(cauchyerrorpost, proposal, start, 50000, darwin$difference)
proposal2 <- list(var = fitlaplace$var, mu = t(fitlaplace$mode))
fitindep <- indepmetrop(cauchyerrorpost, proposal2, start, 50000, darwin$difference)
fitgibbs <- gibbs(cauchyerrorpost, start, 50000, c(12, 0.75), darwin$difference)

apply(fitrw$par, 2, mean)

apply(fitrw$par, 2, sd)

# Table 6.2. Summaries of the marginal posterior densities of μ and log(σ) using five computational methods.

# 6.10 Analysis of the Stanford Heart Transplant Data

str(stanfordheart)
#data(stanfordheart)

start <- c(0, 3, -1)
laplacefit <- laplace(transplantpost, start, stanfordheart)
laplacefit

proposal <- list(var = laplacefit$var, scale = 2)
s <- rwmetrop(transplantpost, proposal, start, 10000, stanfordheart)
s$accept

par(mfrow = c(2, 2))
tau <- exp(s$par[ , 1])
plot(density(tau), main = "TAU")
lambda <- exp(s$par[ , 2])
plot(density(lambda), main = "LAMBDA")
p <- exp(s$par[ , 3])
plot(density(p), main = "P")

# Fig. 6.11. Posterior densities of parameters τ, λ, and p in the Pareto survival model.

apply(exp(s$par), 2, quantile, c(0.05, 0.5, 0.95))

par(mfrow = c(1, 1))
t <- seq(1, 240)
p5 <- 0*t
p50 <- 0*t
p95 <- 0*t
for (j in 1:240) { 
  S <- (lambda/(lambda + t[j]))^p
  q <- quantile(S, c(0.05, 0.5, 0.95))
  p5[j] <- q[1]
  p50[j] <- q[2]
  p95[j] <- q[3]
}
plot(t, p50, type = "l", ylim = c(0, 1), ylab = "Prob(Survival)", xlab = "time")
lines(t, p5 , lty = 2)
lines(t, p95, lty = 2)

# Fig. 6.12. Posterior distribution of probability of survival S(t) for heart transplant
# patients. Lines correspond to the 5th, 50th, and 95th percentiles of the posterior of
# S(t) for each time t.

# 6.11 Further Reading

# 6.12 Summary of R Functions

# 6.13 Exercises

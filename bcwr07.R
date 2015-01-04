# Bayesian Computation with R (Second Edition)
# by Jim Albert

library(LearnBayes)
library(lattice)

# 7 Hierarchical Modeling

# 7.1 Introduction

# 7.2 Three Examples

str(sluggerdata)

# fit logistic model for home run data for a particular player

logistic.fit <- function(player) {
  d <- subset(sluggerdata, Player == player)
  x <- d$Age
  x2 <- d$Age^2
  response <- cbind(d$HR, d$AB - d$HR)
  list(Age = x, p = glm(response ~ x + x2, family = binomial)$fitted)
}

names <- unique(sluggerdata$Player)
newdata <- NULL
for (j in 1:9) {
  fit <- logistic.fit(as.character(names[j]))
  newdata <- rbind(newdata, data.frame(as.character(names[j]), fit$Age, fit$p))
}
names(newdata) <- c("Player", "Age", "Fitted")
xyplot(Fitted ~ Age | Player, data = newdata, type = "l", lwd = 3, col = "black", layout = c(3, 3))

# Fig. 7.1. Plots of fitted career trajectories for nine baseball players as a function of their age.

# 7.3 Individual and Combined Estimates

str(hearttransplants)

with(hearttransplants, plot(log(e), y/e, xlim = c(6, 9.7), xlab = "log(e)", ylab = "y/e"))
with(hearttransplants, text(log(e), y/e, labels = as.character(y), pos = 4))

# Fig. 7.2. Plot of death rates against log exposure for all hospitals. Each point is
# labeled by the number of observed deaths.

# 7.4 Equal Mortality Rates?

sum(hearttransplants$y)
sum(hearttransplants$e)

lambda <- rgamma(1000, shape = 277, rate = 294681)

ys94 <- rpois(1000, hearttransplants$e[94]*lambda)

hist(ys94, breaks = seq(0.5, max(ys94) + 0.5))
with(hearttransplants, lines(c(y[94], y[94]), c(0, 120), lwd = 3))

# Fig. 7.3. Histogram of simulated draws from the posterior predictive distribution
# of y94*. The actual number of transplant deaths is shown by a vertical line.

lambda <- rgamma(1000, shape = 277, rate = 294681)
prob.out <- function(i) {
  ysi <- rpois(1000, hearttransplants$e[i]*lambda)
  pleft  <- sum(ysi <= hearttransplants$y[i])/1000
  pright <- sum(ysi >= hearttransplants$y[i])/1000
  min(pleft, pright)
}
pout <- sapply(1:94, prob.out)

plot(log(hearttransplants$e), pout, ylab = "Prob(extreme)")

# Fig. 7.4. Scatterplot of predictive probabilities of “at least as extreme” against log
# exposures for all observations.

# 7.5 Modeling a Prior Belief of Exchangeability

pgexchprior <- function(lambda, pars) {
  alpha <- pars[1]
  a <- pars[2]
  b <- pars[3]
  (alpha - 1)*log(prod(lambda)) - (2*alpha + a)*log(alpha*sum(lambda) + b)
}

alpha <- c(5, 20, 80, 400)
par(mfrow = c(2, 2))
for (j in 1:4)
  mycontour(pgexchprior, c(0.001, 5, 0.001, 5), c(alpha[j], 10, 10),
            main = paste("ALPHA = ", alpha[j]), xlab = "LAMBDA 1", ylab = "LAMBDA 2")

# Fig. 7.5. Contour graphs of the exchangeable prior on (λ1 , λ2 ) when μ has an inverse
# gamma(10, 10) distribution and for values of the precision parameter α = 5, 20, 80, and 400.

# 7.6 Posterior Distribution

# 7.7 Simulating from the Posterior

datapar <- list(data = hearttransplants, z0 = 0.53)
start <- c(2, -7)
fit <- laplace(poissgamexch, start, datapar)
fit

par(mfrow = c(1, 1))
mycontour(poissgamexch, c(0, 8, -7.3, -6.6), datapar, xlab = "log alpha", ylab = "log mu")

# Fig. 7.6. Contour plot of the posterior density of (log α, log μ) for the heart transplant
# example. Contour lines are drawn at 10%, 1%, and .1% of the modal value.

start <- c(4, -7)
fitgibbs <- gibbs(poissgamexch, start, 1000, c(1, 0.15), datapar)
fitgibbs$accept

mycontour(poissgamexch, c(0, 8, -7.3, -6.6), datapar, xlab = "log alpha", ylab = "log mu")
points(fitgibbs$par[ , 1], fitgibbs$par[ , 2])

# Fig. 7.7. Contour plot of the posterior density of (log α, log μ) for the heart transplant
# example with a sample of simulated values placed on top.

plot(density(fitgibbs$par[ , 1], bw = 0.2), main = "")
abline(h = 0, col = "grey")

# Fig. 7.8. Density estimate of simulated draws from the marginal posterior of log α.

alpha <- exp(fitgibbs$par[ , 1])
mu <- exp(fitgibbs$par[ , 2])
lam1 <- with(hearttransplants, rgamma(1000, y[1] + alpha, e[1] + alpha/mu))

alpha <- exp(fitgibbs$par[ , 1])
mu <- exp(fitgibbs$par[ , 2])
with(hearttransplants, plot(log(e), y/e, pch = as.character(y)))
for (i in 1:94) {
  lami <- with(hearttransplants, rgamma(1000, y[i] + alpha, e[i] + alpha/mu))
  probint <- quantile(lami, c(0.05, 0.95))
  lines(log(hearttransplants$e[i]) * c(1, 1), probint)
}

# Fig. 7.9. Plot of observed death rates against log exposures together with intervals
# representing 90% posterior probability bands for the true rates {λi}.

# 7.8 Posterior Inferences

# 7.8.1 Shrinkage

shrink <- function(i) mean(alpha/(alpha + hearttransplants$e[i] * mu))
shrinkage <- sapply(1:94, shrink)
plot(log(hearttransplants$e), shrinkage, ylim = c(0, 1))

# Fig. 7.10. Plot of the posterior shrinkages against the log exposures for the heart transplant example.

# 7.8.2 Comparing Hospitals

mrate <- function(i) 
  with(hearttransplants, mean(rgamma(1000, y[i] + alpha, e[i] + alpha/mu)))
hospital <- 1:94
meanrate <- sapply(hospital, mrate)
hospital[meanrate == min(meanrate)]

sim.lambda <- function(i) 
  with(hearttransplants, rgamma(1000, y[i] + alpha, e[i] + alpha/mu))
LAM <- sapply(hospital, sim.lambda)

compare.rates <- function(x) {
  nc <- NCOL(x)
  ij <- as.matrix(expand.grid(1:nc, 1:nc))
  m <- as.matrix(x[ , ij[ , 1]] > x[ , ij[ ,2]]) 
  matrix(colMeans(m), nc, nc, byrow = TRUE)
}
better <- compare.rates(LAM)

better[1:24, 85]

# 7.9 Bayesian Sensitivity Analysis

sir.old.new <- function(theta, prior, prior.new) {
  log.g <- log(prior(theta))
  log.g.new <- log(prior.new(theta))
  wt <- exp(log.g.new - log.g - max(log.g.new - log.g))
  probs <- wt/sum(wt)
  n <- length(probs)
  indices <- sample(1:n, size = n, prob = probs, replace = TRUE)
  theta[indices]
}

prior <- function(theta)
  0.53*exp(theta)/(exp(theta) + 0.53)^2
prior.new <- function(theta)
  5*exp(theta)/(exp(theta) + 5)^2

log.alpha <- fitgibbs$par[ , 1]
log.alpha.new <- sir.old.new(log.alpha, prior, prior.new)

# drawing figure

draw.graph <- function() {
  LOG.ALPHA <- data.frame("prior", log.alpha)
  names(LOG.ALPHA) <- c("Prior", "log.alpha")
  LOG.ALPHA.NEW <- data.frame("new.prior", log.alpha.new)
  names(LOG.ALPHA.NEW) <- c("Prior", "log.alpha")
  D <- densityplot(~log.alpha, group = Prior, 
                   data = rbind(LOG.ALPHA, LOG.ALPHA.NEW), 
                   plot.points = FALSE, 
                   main = "Original Prior and Posterior (solid), New Prior and Posterior (dashed)",
                   lwd = 4, adjust = 2, lty = c(1, 2), 
                   xlab = "log alpha", xlim = c(-3, 5), col = "black")
  update(D, panel = function(...) {
    panel.curve(prior(x)    , lty = 1, lwd = 2, col = "black")
    panel.curve(prior.new(x), lty = 2, lwd = 2, col = "black")
    panel.densityplot(...)
  })
}

draw.graph()

# Fig. 7.11. Prior and posterior distributions of log α for the prior parameter choices
# z0 = 0.53 and z0 = 5 for the heart transplant problem.

# does posterior on true rates depend on prior?

prob.int.rate <- function(i, log.alpha) {
  lami <- with(hearttransplants, rgamma(1000, y[i] + exp(log.alpha), e[i] + exp(log.alpha)/mu))
  quantile(lami, c(0.05, 0.95))
}

ind <- c(1, 25, 57, 89)

ci <- sapply(ind, prob.int.rate, log.alpha)
matplot(matrix(log(hearttransplants$e[ind]), 2, length(ind), byrow = TRUE)       , ci, 
        type = "l", lwd = 4, lty = 1, 
        xlim = c(6, 9.7), ylim = c(0, 0.002), ylab = "True Rate", xlab = "log(e)")

ci <- sapply(ind, prob.int.rate, log.alpha.new)
matplot(matrix(log(hearttransplants$e[ind]), 2, length(ind), byrow = TRUE) + 0.05, ci, 
        type = "l", lwd = 4, lty = 3, add = TRUE)

# 7.10 Posterior Predictive Model Checking

lam94 <- with(hearttransplants, rgamma(1000, y[94] + alpha, e[94] + alpha/mu))

ys94 <- rpois(1000, hearttransplants$e[94]*lam94)

hist(ys94, breaks = seq(-0.5, max(ys94) + 0.5))
lines(hearttransplants$y[94]*c(1, 1), c(0, 100), lwd = 3)

# Fig. 7.12. Histogram of the posterior predictive distribution of y94* for hospital 94
# from the exchangeable model. The observed value of y94 is indicated by the vertical line.

prob.out <- function(i) {
  lami <- with(hearttransplants, rgamma(1000, y[i] + alpha, e[i] + alpha/mu))
  ysi <- rpois(1000, hearttransplants$e[i]*lami)
  pleft  <- sum(ysi <= hearttransplants$y[i])/1000
  pright <- sum(ysi >= hearttransplants$y[i])/1000
  min(pleft, pright)
}
pout.exchange <- sapply(hospital, prob.out)

par(pty = "s")
plot(pout, pout.exchange, xlab = "P(extreme), equal means", ylab = "P(extreme), exchangeable")
abline(0, 1)

# Fig. 7.13. Scatterplot of posterior predictive probabilities of “at least as extreme”
# for the equal means and exchangeable models.

# 7.11 Further Reading

# 7.12 Summary of R Functions

# 7.13 Exercises

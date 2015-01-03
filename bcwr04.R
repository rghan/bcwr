# Bayesian Computation with R (Second Edition)
# by Jim Albert

library(LearnBayes)

# 4 Multiparameter Models

# 4.1 Introduction

# 4.2 Normal Data with Both Parameters Unknown

str(marathontimes)
d <- mycontour(normchi2post, c(220, 330, 500, 9000), marathontimes$time, xlab = "mean", ylab = "variance")

S <- with(marathontimes, sum((time - mean(time))^2))
n <- length(marathontimes$time)
sigma2 <- S/rchisq(1000, n - 1)
mu <- rnorm(1000, mean = mean(marathontimes$time), sd = sqrt(sigma2)/sqrt(n))

points(mu, sigma2)

# Fig. 4.1. Contour plot of the joint posterior distribution of (μ,σ2) for a normal 
# sampling model. The points represent a simulated random sample from this distribution.

quantile(mu, c(0.025, 0.975))

quantile(sqrt(sigma2), c(0.025, 0.975))


# 4.3 A Multinomial Model

alpha <- c(728, 584, 138)
theta <- rdirichlet(1000, alpha)

hist(theta[, 1] - theta[, 2], main = "")

# Fig. 4.2. Histogram of simulated sample of the marginal posterior distribution of
# θ1 − θ2 for the multinomial sampling example.

str(election.2008)
#data(election.2008)
#attach(election.2008)

prob.Obama <- function(j) {
  p <- with( election.2008, rdirichlet(5000, 500*c(M.pct[j], O.pct[j], 100 - M.pct[j] - O.pct[j])/100 + 1))
  mean(p[ , 2] > p[ , 1])
}

Obama.win.probs <- sapply(1:51, prob.Obama)

sim.election <- function() {
  winner <- rbinom(51, 1, Obama.win.probs)  
  sum(election.2008$EV*winner)         
}

sim.EV <- replicate(1000, sim.election())

hist(sim.EV, min(sim.EV):max(sim.EV), col = "blue")
abline(v = 365, lwd = 3)  # Obama received 365 votes
text(375, 30, "Actual \n Obama \n total")

# Fig. 4.3. Histogram of 1000 simulated draws of the total electoral vote for Barack
# Obama in the 2008 U.S. presidential election. The actual electoral vote of 365 is
# indicated by a vertical line.

# 4.4 A Bioassay Experiment

x <- c(-0.86, -0.3, -0.05, 0.73)
n <- c(5, 5, 5, 5)
y <- c(0, 1, 3, 5)

data.frame(Dose = x, Deaths = y, SampleSize = n)

# Table 4.1. Data from the bioassay experiment.

data <- cbind(x, n, y)

glmdata <- cbind(y, n - y)
results <- glm(glmdata ~ x, family = binomial)
summary(results)

# when x = -0.7, median and 90th percentile of p are (0.2, 0.40)
# when x = +0.6, median and 90th percentile of p are (0.8, 0.95)

a1.b1 <- beta.select(list(p = 0.5, x = 0.2), list(p = 0.9, x = 0.50))
a1.b1
a2.b2 <- beta.select(list(p = 0.5, x = 0.8), list(p = 0.9, x = 0.98))
a2.b2

data.frame(Dose = c(-0.7, 0.6), Deaths = c(1.12, 2.10), SampleSize = c(4.68, 2.84))

# Table 4.2. Prior information in the bioassay experiment.

prior <- rbind(c(-0.7, 4.68, 1.12), c(0.6, 2.10, 0.74))
data.new <- rbind(data, prior)

# plot prior

plot(c(-1, 1), c(0, 1), type = "n", xlab = "Dose", ylab = "Prob(death)")
lines(-0.7*c(1, 1), qbeta(c(0.25, 0.75), a1.b1[1], a1.b1[2]), lwd = 4)
lines( 0.6*c(1, 1), qbeta(c(0.25, 0.75), a2.b2[1], a2.b2[2]), lwd = 4)
points(c(-0.7, 0.6), qbeta(0.5, c(a1.b1[1], a2.b2[1]), c(a1.b1[2], a2.b2[2])), pch = 19, cex = 2)
text(-0.3, 0.2, "Beta(1.12, 3.56)")
text( 0.2, 0.8, "Beta(2.10, 0.74)")
response <- rbind(a1.b1, a2.b2)
x <- c(-0.7, 0.6)
fit <- glm(response ~ x, family = binomial)
curve(exp(fit$coef[1] + fit$coef[2]*x)/(1 + exp(fit$coef[1] + fit$coef[2]*x)), add = T)

# Fig. 4.4. Illustration of conditional means prior for the bioassay example. In each 
# bar, the point corresponds to the median and the endpoints correspond to the quartiles
# of the prior distribution for each beta distribution.

mycontour(logisticpost, c(-3, 3, -1, 9), data.new, xlab = "beta0", ylab = "beta1")

# Fig. 4.5. Contour plot of the posterior distribution of (β0,β1) for the bioassay
# example. The contour lines are drawn at 10%, 1%, and .1% of the model height.

s <- simcontour(logisticpost, c(-2, 3, -1, 11), data.new, 1000)
points(s)

# Fig. 4.6. Contour plot of the posterior distribution of (β0, β1) for the bioassay
# example. A simulated random sample from this distribution is shown on top of the
# contour plot.

plot(density(s$y), xlab = "beta1", main = "")
abline(h = 0, col = "grey")

# Fig. 4.7. Density of simulated values from the posterior of the slope parameter β1 in the bioassay example.

theta <- -s$x/s$y
hist(theta, xlab = "LD-50", breaks = 20)

# Fig. 4.8. Histogram of simulated values of the LD-50 parameter −β0/β1 in the bioassay example.

quantile(theta, c(0.025, 0.975))

# 4.5 Comparing Two Proportions

sigma <- c(2, 1, 0.5, 0.25)
plo <- 0.0001
phi <- 0.9999
par(mfrow = c(2, 2))
for (i in 1:4)
  mycontour(howardprior, c(plo, phi, plo, phi), c(1, 1, 1, 1, sigma[i]),
            main = paste("sigma = ", as.character(sigma[i])),
            xlab = "p1", ylab = "p2")

# Fig. 4.9. Contour graphs of Howard’s dependent prior for values of the association parameter σ = 2, 1, 0.5, and 0.25.

# Table 4.3. Pearson’s example.

sigma <- c(2, 1, 0.5, 0.25)
par(mfrow = c(2, 2))
for (i in 1:4) {
  mycontour(howardprior, c(plo, phi, plo, phi),
            c(1 + 3, 1 + 15, 1 + 7, 1 + 5, sigma[i]),
            main = paste("sigma = ", as.character(sigma[i])),
            xlab = "p1", ylab = "p2")
  lines(c(0, 1), c(0, 1))
}

# Fig. 4.10. Contour graphs of the posterior for Howard’s dependent prior for values
# of the association parameter σ = 2, 1, 0.5, and 0.25.

s <- simcontour(howardprior, c(plo, phi, plo, phi), c(1 + 3, 1 + 15, 1 + 7, 1 + 5, 2), 1000)
sum(s$x > s$y)/1000

# Table 4.4. Posterior probabilities of the hypothesis.

# 4.6 Further Reading

# 4.7 Summary of R Functions

# 4.8 Exercises

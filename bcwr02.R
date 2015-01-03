# Bayesian Computation with R (Second Edition)
# by Jim Albert

# 2 Introduction to Bayesian Thinking

# 2.1 Introduction

# 2.2 Learning About the Proportion of Heavy Sleepers

# 2.3 Using a Discrete Prior 

p <- seq(0.05, 0.95, by = 0.1)
prior <- c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0)
prior <- prior/sum(prior)
plot(p, prior, type = "h", ylab = "Prior Probability")

# Fig. 2.1. A discrete prior distribution for a proportion p.

data <- c(11, 16)
post <- pdisc(p, prior, data)
round(cbind(p, prior, post), 2)

library(lattice)
PRIOR <- data.frame("prior", p, prior)
POST <- data.frame("posterior", p, post)
names(PRIOR) <- c("Type", "P", "Probability")
names(POST) <- c("Type", "P", "Probability")
data <- rbind(PRIOR, POST)
xyplot(Probability ~ P | Type, data = data, layout = c(1, 2), type = "h", lwd = 3, col = "black")

# Fig. 2.2. Prior and posterior distributions for a proportion p using a discrete prior.

# 2.4 Using a Beta Prior 

quantile2 <- list(p = 0.9, x = 0.5)
quantile1 <- list(p = 0.5, x = 0.3)
beta.select(quantile1, quantile2)

a <- 3.26
b <- 7.19
s <- 11
f <- 16
curve(dbeta(x, a + s, b + f), from = 0, to = 1, xlab = "p", ylab = "Density", lty = 1, lwd = 4)
curve(dbeta(x, s + 1, f + 1), add = TRUE, lty = 2, lwd = 4)
curve(dbeta(x, a, b), add = TRUE, lty = 3, lwd = 4)
legend(0.7, 4, c("Prior", "Likelihood", "Posterior"), lty = c(3, 2, 1), lwd = c(3, 3, 3))

# Fig. 2.3. The prior density g(p), the likelihood function L(p), and the posterior 
# density g(p|data) for learning about a proportion p.

1 - pbeta(0.5, a + s, b + f)

qbeta(c(0.05, 0.95), a + s, b + f)

ps <- rbeta(1000, a + s, b + f)

hist(ps, xlab = "p",main = "")

# Fig. 2.4. A simulated sample from the beta posterior distribution of p.

sum(ps >= 0.5)/1000

quantile(ps, c(0.05, 0.95))

# 2.5 Using a Histogram Prior 

midpt <- seq(0.05, 0.95, by = 0.1)
prior <- c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0)
prior <- prior/sum(prior)
curve(histprior(x, midpt, prior), from = 0, to = 1, ylab = "Prior density", ylim = c(0, 0.3))

# Fig. 2.5. A histogram prior for a proportion p.

curve(histprior(x, midpt, prior) * dbeta(x, s + 1, f + 1), from = 0, to = 1, ylab = "Posterior density")

# Fig. 2.6. The posterior density for a proportion using a histogram prior

p <- seq(0, 1, length = 500)
post <- histprior(p, midpt, prior) * dbeta(p, s + 1, f + 1)
post <- post/sum(post)

ps <- sample(p, replace = TRUE, prob = post)

hist(ps, xlab = "p", main = "")

# Fig. 2.7. A histogram of simulated draws from the posterior distribution of p with the use of a histogram prior.

# 2.6 Prediction

p <- seq(0.05, 0.95, by = 0.1)
prior <- c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0)
prior <- prior/sum(prior)
m <- 20
ys <- 0:20
pred <- pdiscp(p, prior, m, ys)
round(cbind(0:20, pred), 3)

ab <- c(3.26, 7.19)
m <- 20
ys <- 0:20
pred <- pbetap(ab, m, ys)

p <- rbeta(1000, 3.26, 7.19)

y <- rbinom(1000, 20, p)

table(y)

freq <- table(y)
ys <- as.integer(names(freq))
predprob <- freq/sum(freq)
plot(ys, predprob, type = "h", xlab = "y", ylab = "Predictive Probability")

# Fig. 2.8. A graph of the predictive probabilities of the number of sleepers y ̃ in a future sample of size 20 when
# the proportion is assigned a beta(3.26, 7.19) prior.

dist <- cbind(ys, predprob)
dist

covprob <- 0.9
discint(dist, covprob)

# 2.7 Further Reading

# 2.8 Summary of R Functions

# 2.9 Exercises

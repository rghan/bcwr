# Bayesian Computation with R (Second Edition)
# by Jim Albert

# 1 An Introduction to R

# 1.1 Overview

# 1.2 Exploring a Student Dataset 

# 1.2.1 Introduction to the Dataset

# 1.2.2 Reading the Data into R 

# studentdata <- read.table("studentdata.txt", sep = "\t", header = TRUE)

# This dataset is also available as part of the LearnBayes package. 

library(LearnBayes)
data(studentdata)

studentdata[1, ]

# attach(studentdata)

# 1.2.3 R Commands to Summarize and Graph a Single Batch 

table(studentdata$Drink)


barplot(table(studentdata$Drink), xlab = "Drink", ylab = "Count")

# Fig. 1.1. Barplot of the drinking preference of the statistics students.

studentdata$hours.of.sleep <- with(studentdata, WakeUp - ToSleep)

summary(studentdata$hours.of.sleep)

hist(studentdata$hours.of.sleep, main = "")

# Fig. 1.2. Histogram of the hours of sleep of the statistics students.

# 1.2.4 R Commands to Compare Batches 

boxplot(hours.of.sleep ~ Gender, data = studentdata, ylab = "Hours of Sleep")

# Fig. 1.3. Parallel boxplots of the hours of sleep of the male and female students.

str(studentdata)

female.Haircut <- with(studentdata, Haircut[Gender == "female"])
male.Haircut <- with(studentdata, Haircut[Gender == "male"])
summary(female.Haircut)
summary(male.Haircut)

# 1.2.5 R Commands for Studying Relationships 


par(pty = "s")

with(studentdata, plot(jitter(ToSleep), jitter(hours.of.sleep)))

# Fig. 1.4. Scatterplot of wake-up time and hours of sleep for students. 

with(studentdata, plot(jitter(ToSleep), jitter(hours.of.sleep)))
fit <- lm(hours.of.sleep ~ ToSleep, data = studentdata)
fit
abline(fit)

# Fig. 1.5. Scatterplot of wake-up time and hours of sleep for students with least- squares line plotted on top.

# 1.3 Exploring the Robustness of the t Statistic 

# 1.3.1 Introduction 

# 1.3.2 Writing a Function to Compute the t Statistic 

x <- rnorm(10, mean = 50, sd = 10)
y <- rnorm(10, mean = 50, sd = 10)

m <- length(x)
n <- length(y)

sp <- sqrt(((m - 1)*sd(x)^2+ (n - 1)*sd(y)^2)/(m + n - 2))

t.stat <- (mean(x) - mean(y))/(sp*sqrt(1/m+ 1/n))
t.stat

tstatistic <- function(x, y) {
  m <- length(x)
  n <- length(y)
  sp <- sqrt(((m - 1)*sd(x)^2 + (n - 1)*sd(y)^2)/(m + n - 2))
  t.stat <- (mean(x) - mean(y))/(sp*sqrt(1/m + 1/n))
  return(t.stat)
}

#source("tstatistic.R")

data.x <- c(1, 4, 3, 6, 5)
data.y <- c(5, 4, 7, 6, 10)
tstatistic(data.x, data.y)

# 1.3.3 Programming a Monte Carlo Simulation

alpha <- 0.1  # sets alpha
m <- 10       # sets m
n <- 10       # sets n
N <- 10000    # sets the number of simulations
n.reject <- 0 # counter of num. of rejections
for (i in 1:N) {
  x <- rnorm(m, mean = 0, sd = 1) # simulates xs from population 1
  y <- rnorm(n, mean = 0, sd = 1) # simulates ys from population 2
  t.stat <- tstatistic(x, y)      # computes the t statistic
  if (abs(t.stat) > qt(1 - alpha/2, n + m - 2))
    n.reject <- n.reject + 1 # reject if |T| exceeds critical pt 
}
true.sig.level <- n.reject/N # est. is proportion of rejections
true.sig.level

# 1.3.4 The Behavior of the True Significance Level Under Different Assumptions 

m <- 10
n <- 10
my.tsimulation <- function() {
  tstatistic(rnorm(m, mean = 10, sd = 2), rexp(n, rate = 1/10))
}

tstat.vector <- replicate(10000, my.tsimulation())

plot(density(tstat.vector), xlim = c(-5, 8), ylim = c(0, 0.4), lwd = 3)
curve(dt(x, df = 18), add = TRUE)
legend(4, 0.3, c("exact", "t(18)"), lwd = c(3,1))

# Fig. 1.6. Exact sampling density of the t statistic when sampling from normal and exponential distributions. 
# The t sampling density assuming normal populations is also displayed.

# 1.4 Further Reading

# 1.5 Summary of R Functions 

# 1.6 Exercises

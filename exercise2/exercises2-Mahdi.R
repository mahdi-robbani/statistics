setwd('exercise2/')
cells <- read.csv("cell_types.csv", na.strings = "")

#Brain cell dataset
##Ex 4.1
left <- cells$ef__peak_t_ramp[cells$specimen__hemisphere == "left"]
left <- left[!is.na(left)]
hist(left)

right <-
  cells$ef__peak_t_ramp[cells$specimen__hemisphere == "right"]
right <- right[!is.na(right)]
hist(right)

##Ex 4.2
plot(density(right), col = "blue", main = "Kernel Density Estimations for left and right hemispheres")
lines(density(left), col = "red")
legend(
  "topright",
  legend = c("Left Hemisphere", "Right Hemisphere"),
  col = c("red", "blue"),
  lty = 1
)
###The right hemisphere has a higher peak

##Ex 4.3
male <-
  cells$ef__peak_t_ramp[cells$donor__species == "Homo Sapiens" &
                          cells$donor__sex == "Male"]
female <-
  cells$ef__peak_t_ramp[cells$donor__species == "Homo Sapiens" &
                          cells$donor__sex == "Female"]
female <- female[!is.na(female)]

plot(density(female), col = "red", main = "Kernel Density Estimations for males and females")
lines(density(male), col = "blue")
legend(
  "topright",
  legend = c("Female", "Male"),
  col = c("red", "blue"),
  lty = 1
)
#Females have a higher peak

#Ex 5
##Ex 5.1
human <- cells$ef__peak_t_ramp[cells$donor__species == "Homo Sapiens"]
human <- human[!is.na(human)]
theoretical_q <- qlnorm(ppoints(human), sdlog = 0.6, meanlog = 2)

plot(theoretical_q, sort(human), 
     xlab = "Theoretical quantiles",
     ylab = "Empirical quantiles")
abline(0,1, col ="red")

##Ex 5.2
spike <- cells$ef__peak_t_ramp
spike <- spike[!is.na(spike)]
hist(spike, prob = TRUE)
curve(dnorm(x, mean = 6.41, sd = 4.32), add = TRUE, col = "red")

##Ex 5.3
theoretical_spike_q <- qnorm(ppoints(spike), sd = 4.32, mean = 6.41)
plot(theoretical_spike_q, sort(spike), 
     xlab = "Theoretical quantiles",
     ylab = "Empirical quantiles")
abline(0,1, col ="red")


#Emperical mean and variance
#Ex 6.1
trueMean <- 100*0.3
trueVar <- 100*0.3*(1-0.3)
trueSD <- sqrt(trueVar)

eMean <- c()
eVar <- c()
eSD <- c()
loop <- c(10, 100, 1000, 10000)

for (i in 1:length(loop)) {
  eMean[i] <- mean(rbinom(loop[i], size = 100, prob = 0.3))
  eVar[i] <- var(rbinom(loop[i], size = 100, prob = 0.3))
  eSD[i] <- sd(rbinom(loop[i], size = 100, prob = 0.3))
}

plot(
  loop,
  eMean,
  col = "blue",
  pch = 1,
  ylim = (c(0, 40)),
  ylab = "Value",
  main = "Mean, Variance and SD of a binomial function"
)
points(loop, eVar, col = "red", pch = 2)
points(loop, eSD, col = "green", pch = 3)
abline(h = trueMean, col = "blue")
abline(h = trueVar, col = "red")
abline(h = trueSD, col = "green")
legend(
  "top",
  legend = c("Mean", "Variance", "SD"),
  col = c("blue", "red", "green"),
  pch = c(1, 2, 3)
)

#Ex 6.2
eDataN <- function(n){
  replicate(1000, rbinom(n, size = 100, prob = 0.3)))
}
eMeanN <- function(n) {
  return(replicate(1000, mean(rbinom(n, size = 100, prob = 0.3))))
}
eVarN <- function(n) {
  return(replicate(1000, var(rbinom(n, size = 100, prob = 0.3))))
}
eSDN <- function(n) {
  return(replicate(1000, sd(rbinom(n, size = 100, prob = 0.3))))
}

eData10 <- eDataN(10)
eMean10 <- eMeanN(10)
eVar10 <- eVarN(10)
eSD10 <- eSDN(10)
hist(eMean10)

#Ex 6.3
standardErrorMean <- sd(rbinom(10, size = 100, prob = 0.3))/sqrt(10-1) #divided by 10 since i used n= 10
sd(eMean10)

#Standard error of the mean <- theoretial value
#standard deviation of the empirical mean <- emperical value
#These two should be very close, ideally the same


#Ex 6.4
hist(eMean10, prob= TRUE, ylim=c(0, 0.3))
curve(dnorm(x, mean = 100*0.3, sd = standardErrorMean), add = TRUE, col="red")
legend("right", col = ("red"), lty = 1, legend = c("True density")  )
  
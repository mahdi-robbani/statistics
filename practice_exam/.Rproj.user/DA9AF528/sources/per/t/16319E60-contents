setwd('exercise3/')
angles <- read.table("angles.txt")
angles <- as.vector(angles$x)

#Ex 1.2
dsink <- function(x, k = 1, lg = FALSE) {
  #returns a vector of sinx^k or log(sin^k)
  sinintegral <- integrate(function(x) sin(x) ^ k, lower = 0, upper = pi)$value
  if (lg == FALSE) {
    return(sin(x)^k/sinintegral)
  }
  else{
    return(log(sin(x)^k/sinintegral))
  }
}

curve(dsink(x, k =10), from =0, to=pi)
integrate(function(x) sin(x) ^ 10, lower = 0, upper = pi)$value
integrate(function(x) dsink(x), lower=0, upper=pi)$value

dsinkMLL <- function(k, xvals) {
  return(-sum(dsink(xvals, k, lg=TRUE)))
}

#Ex 1.3
optimize(f = dsinkMLL, xvals = angles, interval = c(0, 100))

#Ex 1.4
hist(
  angles,
  xlim = c(0.5, 2.5),
  probability = TRUE
)
mink <- optimize(f = dsinkMLL, xvals = angles, interval = c(0, 12))$minimum
curve(dsink(x, k = mink), add = TRUE, col = "blue")

#Ex 2
##Ex 2.1
isidata <- read.table("neuronspikes.txt")
isidata <- isidata$V1
ISImean <- mean(isidata)
lambda <- 1/ISImean
lambda

##Ex 2.2
gammaMLL <- function(parameters, xvals){
  return(-sum(dgamma(xvals, shape = parameters[1], rate = parameters[2], log = TRUE)))
}

parametersGamma <- optim(par = c(1, 1), fn = gammaMLL, xvals = isidata)$par
paste("alpha =", parametersGamma[1], "beta =", parametersGamma[2])

#Ex 2.3
beta <- mean(isidata)/var(isidata)
beta

#Ex 3.4
dInvNorm <- function(x, mu, lambda, lg = FALSE){
  if(lg == TRUE){
    return(log(sqrt(lambda/(2*pi*(x^3)))*exp(-lambda*((x-mu)^2)/(2*(mu^2)*x))))
  }
  else{
    return(sqrt(lambda/(2*pi*(x^3)))*exp(-lambda*((x-mu)^2)/(2*(mu^2)*x)))
  }
}

invNormMLH <- function(parameters, xvals){
  return(-sum(dInvNorm(xvals, mu=parameters[1], lambda = parameters[2], lg = TRUE)))
}

parametersinvNorm <- optim(par = c(0.5, 0.5), fn = invNormMLH, xvals = isidata)$par
paste("Mu = ", parametersinvNorm[1])
paste("Lambda =", parametersinvNorm[2])

hist(
  isidata,
  xlim = c(0, 6),
  ylim = c(0, 1.2),
  probability = TRUE
)
curve(dInvNorm(x, parametersinvNorm[1], parametersinvNorm[2]), col = "blue", add = TRUE)

#ex4.1
cells <- read.csv("cell_types.csv", na.strings = "")
rampSpike <- cells$ef__peak_t_ramp
rampSpike <- rampSpike[!is.na(rampSpike)]

dlnormMLL <- function(parameters, xvals){
  return(-sum(dlnorm(xvals, meanlog = parameters[1], sdlog = parameters[2], log = TRUE)))
}

parametersDlnorm <- optim(par = c(1, 1), fn = dlnormMLL, xvals = rampSpike)$par

paste("Meanlog = ", parametersDlnorm[1])
paste("SDlog =", parametersDlnorm[2])

#x4.2
logRampSpike <- log(rampSpike)

dnormMLL <- function(parameters, xvals){
  return(-sum(dnorm(xvals, mean = parameters[1], sd = parameters[2], log = TRUE)))
}

parametersDnorm <- optim(par = c(1, 1), fn = dnormMLL, xvals = logRampSpike)$par

?dnorm
paste("Mean = ", parametersDnorm[1])
paste("SD =", parametersDnorm[2])

#Ex 4.3
rampSpike <- cells$ef__peak_t_ramp
maleRampSpike <- rampSpike[cells$donor__species == "Homo Sapiens" & cells$donor__sex == "Male"]
femaleRampSpike <- rampSpike[cells$donor__species == "Homo Sapiens" & cells$donor__sex == "Female"]
femaleRampSpike <- femaleRampSpike[!is.na(femaleRampSpike)]

parametersDlnormMale <- optim(par = c(1, 1), fn = dlnormMLL, xvals = maleRampSpike)$par
parametersDlnormFemale <- optim(par = c(1, 1), fn = dlnormMLL, xvals = femaleRampSpike)$par

curve(dInvNorm(x, parametersDlnormMale[1], parametersDlnormMale[2]), col = "blue", ylab = "Density")
curve(dInvNorm(x, parametersDlnormFemale[1], parametersDlnormFemale[2]), col = "red", add = TRUE)
legend("topright",legend = c("Male", "Female"), col = c("blue", "red"), lty = c(1, 1))

#Ex 5.1
JKprob <- function(X, Y, alpha, t=1){
  if(X==Y){
    return(0.25+0.75*exp(-4*alpha*t))
  }
  else{
    return(0.25-0.25*exp(-4*alpha*t))
  }
}

#5.3

#5.4
#Implement probability function
distJK <- function(pair, alpha){
  A <- sapply(JK_pairs,function(x) x[1]) #select first element of each vector inside list
  B <- sapply(JK_pairs,function(x) x[2])#select second element (list of letters)
  n <- length(pair)
  n1 <- sum(A == B)
  n2 <- n - n1
  match <- n1*log(JKprob(X="X", Y="X", alpha, t)*0.25)
  mismatch <- n2*log(JKprob(X="X", Y="Y", alpha, t)*0.25)
  return(match+mismatch)
}

#Implement sampling procedure
simulate_JKpair <- function(start="A", alpha=2, t = 10){
  mutations <- rpois(1, lambda = 4*alpha*t)
  if(mutations == 0){
    end <- start
    return(c(start, end))
  } 
  else {
    for (i in 1:mutations){
      end <- sample(c("A", "C", "T", "G"), size = 1)
    }
  }
  return(c(start, end))
}

JK_pairs <- list()
for (i in 1:1000){
  JK_pairs[[i]] <- simulate_JKpair(t=1000, alpha = 1e-5) 
}



#Minus log likelihood function
distJKMLL <- function(alpha, xvals, t = 1){
  probs <- sapply(xvals, function(x) JKprob(X = x[1], Y = x[2], alpha = alpha, t = t))
  return(-sum(log(probs)))
}

#Solving numerically
optimize(f = distJKMLL, xvals = JK_pairs, interval = c(0, 0.01), t = 1)


####Gherardo solution
JK_table <- data.frame(x = sapply(JK_pairs,function(x) x[1]), y= sapply(JK_pairs,function(x) x[2]))

f_jk <- function(x, y, t, a){
  if (x == y){
    0.25 + 0.75 * exp(-4*a*t)
  }else{
    0.25 - 0.25 * exp(-4*a*t)
  }
}


mll_jk <- function(a, data, t = 1){
  -sum(sapply(1:nrow(data), function(i){
    log(f_jk(x = data$x[i], y = data$y[i], t = t, a = a))
  }))
}

mll_jk(1, data = JK_table)

a_est <- optimize(f = mll_jk, interval = c(0,10), data = JK_table, t = 1000)$minimum

a_est




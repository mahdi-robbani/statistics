load("exam1819.RData")

#Question 1.1
probability_mice <- as.data.frame(table(radiation))
probability_mice$prob <- probability_mice$Freq/nrow(radiation)
probability_mice[probability_mice$dead == 0,]

death <- probability_mice[probability_mice$dead == 0,]
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
death$dose <- as.numeric.factor(death$dose) #The correct approach is to divide the contignency tables, e.g dividing mice with 201 dose and no treatmeant that died, with total number of 201 dose mice that receivieced no treatment. You divided by all mice

strep <- death[death$treatment == 0,]
strep <- strep[,c("dose", "prob")]

cont <- death[death$treatment == 1,]
cont <- cont[,c("dose", "prob")]

plot(strep$dose, strep$prob, type = "o", col = "red", xlab = "Dose", ylab = "Probability", main = "Probability of Death", ylim= c(0,0.2))
lines(cont$dose, cont$prob, type = "o", col = "blue")
legend("topright", c("Streptomycin", "Control"), lty = 1, col = c("red", "blue"))

########Gherardo solution
table(radiation)[,,2] # select second contigency table, commas dictate dimensions, i.e. 3 values separated by commas for a 3 dimensional contingency table

probs_death <- table(radiation)[,,2] /  
  (table(radiation)[,,1] + table(radiation)[,,2])  
probs_death

doses <- dimnames(probs_death)$dose
plot(doses, probs_death[, 1], type = "b", col = "red", ylim = c(0,1), 
     ylab = "probability", main = "Probability of death")
lines(doses, probs_death[, 2], type = "b", col = "blue")
legend("topleft", legend = c("Streptomycin", "Saline control"), 
       col = c("blue", "red"), lty = 1, pch = 1)

###############

#Question 1.2
chisq.test(table(radiation[radiation$dose == 201, c(2, 3)]))
#p value greater than 0.05 so we cannot reject the null hypothesis. So we cannot reject that death and radiation are indepedant
chisq.test(table(radiation[radiation$dose == 220, c(2, 3)]))
#p value is very small so we reject the null hypothesis. Death is dependant on radiation
chisq.test(table(radiation[radiation$dose == 243, c(2, 3)]))
#p value is very small so we reject the null hypothesis. Death is dependant on radiation
chisq.test(table(radiation[radiation$dose == 260, c(2, 3)]))
#p value is 0.4 so we reject the null hypothesis. Death is dependant on radiation

#Question 1.3
#######################
Wald.stats <- t(sapply(doses, function(d){
  probs <- probs_death[d,]
  n1 <- sum(radiation$dose == d & radiation$treatment == 1)
  n0 <- sum(radiation$dose == d & radiation$treatment == 1) ## in this case they are equal
  delta <- probs[2] - probs[1] ## the order is inverted
  se <- sqrt(probs[1]*(1 - probs[1]) / n0 + probs[2] * (1 - probs[2]) / n1)
  return(list(delta = delta, se = se, W = delta / se))  
}))
Wald.stats

sapply(doses, function(d) {
  probs <- probs_death[d,]
  probs[2]
  })

alpha <- 0.05
z <- qnorm(alpha)
t(sapply(doses, function(d){
  W <- Wald.stats[d,]$W
  return(list(W = W, reject = W < z, pvalue = pnorm(W)))
}))


####################

#Question 1.4
head(dead_dose_logit)
dead_dose_logit <- radiation[radiation$dead == 0,]
dead_dose_logit <- dead_dose_logit[dead_dose_logit$treatment == 1, ]
dead_dose_logit$dead <- NULL

head(radiation[radiation$treatment == 1,])
strep

log_regr_1 <- glm(prob ~ dose, data = strep, family = "binomial")
summary(strep)
strep
plot(predict(log_regr_1))

##############Gherardo
logregr1 <- glm(dead ~ dose, family = binomial, data = radiation[radiation$treatment == 1,])
summary(logregr1)


doses <- as.numeric(doses) ##we need numeric values now
predicted_probs <- predict(logregr1, type = "response", newdata = data.frame(dose = doses)) #type = response inverts the family = binomial (inverse of the link function)

### now we just repeta the plot of Q1.1
plot(doses, probs_death[, 2], type = "b", col = "blue", ylim = c(0,1), 
     ylab = "probability", main = "Probability of death (Streptomycin)")
lines(doses, predicted_probs, type = "b", col = "purple")
legend("topleft", legend = c("Empirical freq.", "Predicted (logregr1)"),
       col = c("blue", "purple"), lty = 1, pch = 1)

#############

#Question 2.1
fat_model <- lm(fat ~ ., data = bodyfat)
summary(fat_model)
#Abdomen circumference is the most relevant followed by wrist and forearm cicumference. This is because these measurements have a small p value, specially abdomen.
#The knee measurement has a pvalue of 0.97930 which is much higher than the significance level of 0.05, so we cannot reject the null hypothesis that the coefficient of knee is 0.

#Question 2.2
fat_model_start <- lm(fat ~ 1, data=bodyfat)
fat_model_forward <- step(fat_model_start, scope= formula(fat_model), direction = "forward")

fat_model_backward <- step(fat_model, direction ="backward")

#Question 2.3
bmi <- bodyfat$weight/(bodyfat$height)^2
bodyfat$bmi <- bmi
bmi_model <- lm(fat ~ bmi + age, data = bodyfat)
summary(bmi_model)
#B_0 is -27.7, B_1 is 1.59 and V_2 is 0.14

#Question 2.4
bootstrap_coefficients <- replicate(1000, expr = {
  bodyfat_bs <- bodyfat[sample(nrow(bodyfat), size = nrow(bodyfat), replace = TRUE),] #sample dataframe
  bmi_bs <- bodyfat_bs$weight/(bodyfat_bs$height)^2
  lm(fat ~ bmi_bs + age, data = bodyfat_bs)$coefficients
})
bootstrap_coefficients <- t(bootstrap_coefficients) #transpose results
bootstrap_coefficients <- as.data.frame(bootstrap_coefficients)
standarderror_bs <- sapply(bootstrap_coefficients, function(x) sd(x))



z <- qnorm(1-0.05/2)
a <- bmi_model$coefficients - z * standarderror_bs
b <- bmi_model$coefficients + z * standarderror_bs
ci_bmi <- c(a, b)
ci_bmi

###Gherardo
x <- t(apply(bootstrap_coefficients, MARGIN = 1, quantile, probs = c(alpha/2, 1 - alpha/2)))
x
confint(bmi_model)

?quantile

#the builtin function had a much smaller interval for both bmi and age.

#Question 2.5
#Haven't done cross validation yet
######### Gherardo
squaredres_cv <- sapply(1:nrow(bodyfat), function(i){
  tmp_step <- lm(formula(fat_model_forward), data = bodyfat[-i,])
  tmp_bmi <- lm(formula(bmi_model), data = bodyfat[-i,])
  prd_step <- predict(tmp_step, newdata = bodyfat[i,])
  prd_bmi <- predict(tmp_bmi, newdata = bodyfat[i,])
  return(c(step = (prd_step - bodyfat$fat[i])^2, bmi =  (prd_bmi - bodyfat$fat[i])^2))
} )
rowSums(squaredres_cv)

#Question 3.1
dgumbel <- function(x, mu=0, b=1){
  z <- (x - mu)/b
  return(1/b * exp(-(z+exp(-z))))
}

pgumbel <- function(q, mu=0, b=1){
  return(exp(-exp(-(q-mu)/b)))
}

qgumbel <- function(p, mu=0, b=1){
  return(mu - b*log(-log(p)))
}

rgumbel <- function(n, mu=0, b=1){
  return(sapply(runif(n), qgumbel, mu = mu, b = b))
}

###test
first <- FALSE
for (mu in -5:5){
  curve(dgumbel(x, mu = mu, b =1), from = -10, to = 10, col = (mu + 5), add = first,
        main = "dgumbel(x, mu, 1)", ylab = "density")
  first <- TRUE
}
legend("right", legend = -5:5, col = -5:5 + 5, lty = 1)

first <- FALSE
for (b in 1:5){
  curve(qgumbel(x, mu = 0, b = b / 5), from = 0, to = 1, col = (b), add = first, 
        main = "qgumbel(x, 0, b)", ylab = "quantile", xlab = "probability", ylim = c(-2, 4))
  first <- TRUE
}
legend("right", legend = (1:5 ) / 5, col = 1:5 + 5, lty = 1)


hist(rgumbel(10000), probability = TRUE, ylim = c(0, 0.4))
curve(dgumbel, add = TRUE, col = "blue")
?hist


#You can check pgumbel and qgumbel by passing a value for pgumbel to get a result and using that value on qgumbel to see if the returned value is the original value. If it did, this willact as a simple sanity check
#Numerical tests include finding the 0th, 50th, 75th and 100th percentile values for the CDF fucntion (qnorm) and then integragrating the dgumbel function to those points and seeing if the appropriate values are obtained, i.e. 0, 0.5, 0.75 and 1.

#Sanity check
qgumbel(0.5)
pgumbel(qgumbel(0.5))

#Numerical check
integrate(dgumbel, lower = -Inf, upper = qgumbel(0.5))



##########Gherardo
res <- (sapply(-100:100, function(mu){
  sapply(1:10, function(b){
    interval <- sort(rnorm(2, mu, 10 * b))
    p <- pgumbel(interval, mu = mu, b = b)
    return(integrate(dgumbel, interval[1], interval[2], mu = mu, b = b)$value - (p[2] - p[1]))
  })
}))
mu <- 10
b <- 2
interval <- sort(rnorm(2, mu, 10 * b))
p <-  pgumbel(interval, mu = mu, b = b)
integrate(dgumbel, interval[1], interval[2], mu = mu, b = b)$value - (p[2] - p[1])


#Question 3.2
b_est <- sqrt(6*var(wind)/pi^2)
mu_est <- mean(wind) - (-digamma(1) * b_est)

hist(wind, probability = TRUE)
curve(dgumbel(x, mu = mu_est, b = b_est), add = TRUE, col = "blue")

qqplot(sort(wind), rgumbel(ppoints(wind), mu= mu_est, b = b_est))
abline(0, 1, col = 'red')

?ppoints

#Question 3.3
gumbelLL <- function(parameters, xvals){
  return(-sum(log(dgumbel(xvals, mu= parameters[1], b = parameters[2]))))
}

parameters_gumbel <- optim(par = c(mu_est, b_est), fn= gumbelLL, xvals= wind)$par

hist(wind, probability = TRUE)
curve(dgumbel(x, mu = mu_est, b = b_est), add = TRUE, col = "blue")
curve(dgumbel(x, mu = parameters_gumbel[1], b = parameters_gumbel[2]), add = TRUE, col = "red")
legend("topright", c("Method of Moments", "Maximum Likelihood"), lty = 1, col = c("blue", "red"))

#Question 3.4
gumbellGaussLL <- function(parameters, xvals){
  return(-sum(dnorm(xvals, mean = parameters[1], sd = parameters[2], log = TRUE)))
}

par_gauss_gumbel <- optim(par = c(0, 1), fn= gumbellGaussLL, xvals= wind)$par

hist(wind, probability = TRUE)
curve(dnorm(x, mean = par_gauss_gumbel[1], sd = par_gauss_gumbel[2]), add = TRUE, col = "blue")

curve(dnorm(x, mean = mean(wind), sd = sd(wind)), add = TRUE, col = "darkgreen") #easer to do this than all the previous steps

#######THIS METHOD IS WRONG
"qqplot(sort(wind), rnorm(ppoints(wind), mean = par_gauss_gumbel[1], sd = par_gauss_gumbel[2]))
abline(0, 1, col = 'red')
#The Gaussian model does not fiy linearly so it is much worse than the Gumbel model"

##USE THIS
qqnorm(wind)
qqline(wind)


###NOT SURE HOW TO DO AIC BIC
##SKIPPED LIKELIHOOD RATIO TEST

#Question 3.5
############Gherado
M <- 1000
pars_bt <- replicate(M, {
  optim(par = c(mu_est, b_est), fn = gumbellGaussLL, data = sample(wind, replace =TRUE))$par
})
SE <- apply(pars_bt, MARGIN = 1, sd)
SE

alpha <- 0.05
z <- qnorm(alpha / 2, lower.tail = FALSE)
cbind("2.5%" = par_mle - SE * z, "97.5%" = par_mle + SE * z)

t(apply(pars_bt, MARGIN = 1, quantile,  probs = c(alpha/2, 1- alpha/2)))


#non parametric bootstrap
wind_para_bs <- replicate(1000, {
  wind_bs <- sample(wind, size = length(wind), replace = TRUE)
  optim(par = c(mu_est, b_est), fn= gumbelLL, xvals= wind_bs)$par
})
wind_para_bs <- t(wind_para_bs)
wind_para_bs <- as.data.frame(wind_para_bs)
head(wind_para_bs)

#standard error
mu_se <- sd(wind_para_bs$V1)
b_se <- sd(wind_para_bs$V2)

#normal quantiles
z <- qnorm(1-0.05/2)
mu_a <- mu_est - z * sqrt(mu_se)
mu_b <- mu_est + z * sqrt(mu_se)
mu_ci <- c(mu_a, mu_b)
mu_ci

z <- qnorm(1-0.05/2)
b_a <- b_est - z * sqrt(b_se)
b_b <- b_est + z * sqrt(b_se)
b_ci <- c(b_a, b_b)
b_ci

#percentile confidence intervals
quantile(wind_para_bs$V1, probs = c(0.05/2, 1 - 0.05/2))
quantile(wind_para_bs$V2, probs = c(0.05/2, 1 - 0.05/2))

#Question 3.6
##Bayesian inference, not done yet, skip
qnorm(1-0.05)

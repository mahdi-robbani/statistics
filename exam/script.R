#1.1
load("exam1920.RData")
hist(obs, breaks = 10, probability = TRUE)
lines(density(obs, 0.5))

plot(ecdf(obs), main = "ECDF plot of obs")
?plot
q <- data.frame(Quantile = c("1st quantile", "1nd quantile", "3rd quantile"),
                Value = quantile(sort(obs))[2:4])
q
###

hist(obs, breaks = 10, probability = TRUE)
for (bw in seq(0.5,2, 0.5)){
  lines(density(obs, bw = bw), col = bw*2  )
}
legend("topright", c("0.5","1", "1.5", "2"), col = seq(0.5,2,0.5)*2, lty = 1)




#1.2
gammaMLH <- function(parameters, xvals){
  return(-sum(dgamma(xvals, shape = parameters[1], rate = parameters[2], log = TRUE)))
}
parametersGamma <- optim(par = c(1, 1), fn = gammaMLH, xvals = obs)$par
valueGamma <- optim(par = c(1, 1), fn = gammaMLH, xvals = obs)$value

hist(obs, breaks = 10, probability = TRUE)
curve(dnorm(x, mean = mean(obs), sd = sd(obs)), add = TRUE, col = "blue")
curve(dlnorm(x, meanlog = mean(log(obs)), sd = sd(log(obs))), add = TRUE, col = "purple")
curve(dgamma(x, shape = parametersGamma[1], rate = parametersGamma[2]), add = TRUE, col = "red")
curve(dexp(x, rate = 1/mean(obs)), add = TRUE, col = "darkgreen")
legend("topright", c("Normal dist", "Log Normal dist", "Gamma dist", "Exp dist"), col = c("blue", "purple", "red", "darkgreen"), lty = 1)

#1.3
qqnorm(obs)
qqline(obs, col = "red")

qqplot(qexp(ppoints(obs), rate = 1/mean(obs)), sort(obs), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", main = "Exponential Q-Q plot")
abline(0, 1, col = "red")

qqnorm(rnorm(100))
qqline(rnorm(100))




loglikeExp <- sum(dexp(obs, rate = 1/mean(obs), log = TRUE))
loglikeGamma <- -valueGamma
delta <- -2 * (loglikeExp - loglikeGamma)
## p-value
pchisq(delta, lower.tail = FALSE, df = 1)

#1.5
models <- list(Normal = dnorm(obs, mean = mean(obs), sd = sd(obs), log= TRUE), Lognormal = dlnorm(obs, meanlog = log(mean(obs)), sd = log(sd(obs)), log = TRUE), Gamma = dgamma(obs, shape = parametersGamma[1], rate = parametersGamma[2], log = TRUE), Exponential = dexp(obs, rate = 1/mean(obs), log = TRUE))

sapply(models, function(x) -2*sum(x) + 2)

z <- qnorm(0.05 / 2, lower.tail = FALSE)
lambda <- 1/(mean(obs))
cbind("2.5%" = lambda - lambda/(sqrt(length(obs))) * z, "97.5%" = lambda + lambda/(sqrt(length(obs))) * z)

se <- lambda/sqrt(length(obs))
wald <- abs((lambda-0.5)/se)
1-pnorm(wald)

wald

#1.7
gamma_parameters_bs <- replicate(1000, {
  obs_bs <- sample(obs, replace = TRUE)
  para_bs <- optim(par = c(1,1), fn = gammaMLH, xvals = obs_bs)$par
  return(c(shape = para_bs[1], rate = para_bs[2]))
})
se <- apply(gamma_parameters_bs, MARGIN = 1, function(x) sd(x))

z <- qnorm(1 - 0.05 / 2)
matrix(parametersGamma + z * se %*% t(c(-1, +1)), dimnames = list(c("shape", "rate"), c("a", "b")), ncol = 2  )

#2.1
head(applejuice)
applejuice_test

applejuice$growth <- as.factor(applejuice$growth)
apple_model <- glm(growth ~ ., family = "binomial", data = applejuice)
summary(apple_model)

update(apple_model, . ~ . + I(ph^2))
#2.3
AIC(apple_model, apple_model2)
BIC(apple_model, apple_model2)
anova(apple_model, apple_model2, test="LRT")


confint(apple_model, level = 0.99)

#2.5
a1_pred <- predict(apple_model, newdata = applejuice_test, type = "response")
a2_pred <- predict(apple_model2, newdata = applejuice_test, type = "response")
data.frame(apple_model = a1_pred, apple_model2 = a2_pred, true_val = applejuice_test$growth_p)

#2.6
new_data <- expand.grid(temp = seq(20:60), nisin = seq(0,80,40))
new_data$ph <- 5
new_data$brix <- 14
grid_pred <- predict(apple_model2, newdata = new_data)
mat_pred <- matrix(grid_pred, nrow=length(seq(20, 60, 20)), byrow = TRUE)
contour(seq(20, 60, 20), seq(0,80,40), mat_pred)
##advanced



data_grid <- expand.grid(temp = 20:60, nisin = 0:80)
data_grid$ph <- 5
data_grid$brix <- 14
pred_grid <- predict(apple_model2, newdata = data_grid, type = "response")
pred_mat <- matrix(pred_grid, nrow=length(20:60), byrow = TRUE)
filled.contour(20:60, 0:80, pred_mat, levels = quantile(pred_mat, probs = seq(0,1,0.05)))


?filled.contour
?contour
#JB
temp <- seq(20, 60, length.out = 100)
nisin <- seq(0, 80, length.out = 100)
data_grid <- expand.grid(temp = temp, nisin = nisin)
data_grid$ph <- 5
data_grid$brix <- 14
pred_grid <- predict(apple_model2, newdata = data_grid, type = "response")
pred_mat <- matrix(pred_grid, nrow=100, byrow = TRUE)
pred_mat[50,]
filled.contour(seq(20, 60, length.out = 100), seq(0, 80, length.out = 100), pred_mat, levels = quantile(pred_mat, probs = seq(0,1,0.1)), xlab= "Temperature", ylab = "Nisin concentration", main = "Contour plot of advanced model")

plot(pred_grid)

dim(pred_mat)
length(20:60)

summary(apple_model2)

## 3.1
hist(brainhead$brainweight)
curve(dnorm(x, mean = mean(brainhead$brainweight), sd = sd(brainhead$brainweight)), add = TRUE, col = "blue")

qqnorm(brainhead$brainweight)
qqline(brainhead$brainweight, col = "red")

##Ex3.2
hist(brainhead$brainweight, probability = TRUE)
curve(dnorm(x, mean = mean(brainhead$brainweight), sd = sd(brainhead$brainweight)), add = TRUE, col = "blue")

#you can do it the analytical way which is better:
z <- qnorm(1 - 0.05 / 2)
a <- mean(brainhead$brainweight) - z * mean(brainhead$brainweight)/sqrt(length(brainhead$brainweight))
b <- mean(brainhead$brainweight) + z * mean(brainhead$brainweight)/sqrt(length(brainhead$brainweight))
conf_mean <- c(a = a, b = b)
conf_mean
mean(brainhead$brainweight)
#You can also do bootstrap but I'm not going to do that since the analytical way is better

#Question 3.3 Consider the following question:
?t.test
head(brainhead)
summary(brainhead)

young <- brainhead$brainweight[brainhead$agerange == 1]
old <- brainhead$brainweight[brainhead$agerange == 2]
t.test(young, old, alternative = "two.sided")

#3.4
young_df <- brainhead[brainhead$agerange == 1,]
old_df <- brainhead[brainhead$agerange == 2,]

brainweight_young <- brainhead$brainweight[brainhead$agerange == 1]
brainweight_old <- brainhead$brainweight[brainhead$agerange == 2]

headsize_young <- brainhead$headsize[brainhead$agerange == 1]
headsize_old <- brainhead$headsize[brainhead$agerange == 2]

head(brainhead)

brain_model <- lm(brainweight ~ headsize, data= brainhead)
summary(brain_model)
brain_model_young <- lm(brainweight ~ headsize, data= young_df)
summary(brain_model_young)
brain_model_old <- lm(brainweight ~ headsize, data= old_df)
summary(brain_model_old)

plot(brainhead$headsize, brainhead$brainweight)
points(headsize_young, brainweight_young, col = "red")
points(headsize_old, brainweight_old, col = "blue")
abline(brain_model)
abline(brain_model_young, col = "red")
abline(brain_model_old, col = "blue")
legend("bottomright", c("Full Data", "Young", "Old"), col = c(1,2, "blue"), lty = 1)

#3.5
qqnorm(residuals(brain_model))
plot(brainhead$headsize, residuals(brain_model))

#3.6
kfold <- function(k, data){
  folds <- list()
  m <- nrow(data) %/% k #group size
  for(i in 1:k){
    if(i == k){
      folds[[i]] <- ((i-1)*m+1):(nrow(data))
    }
    else{
      folds[[i]] <- ((i-1)*m + 1):(i*m)
    }
  }
  return(folds)
}

k <- 10
folds <- kfold(k, brainhead)
folds

leave10 <- function(folds, data){
  output <- c()
  for(i in 1:k){
    fold_i <- sample(folds, replace = FALSE)
    #fit model with all but i row
    model_i <- lm("brainweight ~ headsize", data = brainhead[-fold_i,])
    # predict value for row i
    prediction_i <- predict(model_i, newdata = brainhead[fold_i,])
    #check if prediction is good or not
    se <- (prediction_i - brainhead$brainweight[fold_i])^2
    ##save it
    output[i] <- mean(se)
  }
  return(mean(output))
}

leave10(folds, brainhead)

#3.7
poly_brain <- lm(brainweight ~ headsize + I(headsize^2), data = brainhead)
anova(brain_model, poly_brain, test="LRT")
AIC(brain_model, poly_brain)

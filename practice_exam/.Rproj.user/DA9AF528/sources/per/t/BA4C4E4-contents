setwd("exercise5/")

#Ex 1.1 Plot the observations in a scatter plot. Which one is the true regression function? Plot the true regression function in red in the same plot together with the observations.
GenXY <- function(n=50, a, b, sd=10){
  X <- rnorm(n, mean = 0, sd = sd)
  Y <- sapply(X, function(x) rnorm(1, mean= a*x + b, sd = sd))
  df <- data.frame("X" = X, "Y" = Y)
  return(df)
}

XY <- GenXY(a = 3,b = 12)
plot(XY$X, XY$Y)
abline(a= 3, b = 12, col="red")

#Ex 1.2 Fit now a linear regression using the lm function. Plot the fitted regression line on top of the previous plot and using a different color (e.g. blue)
lm_ab <- lm(Y ~ X, data = XY)
abline(lm_ab, col="blue")

#Ex 1.3 Use the function summary to obtain informations on the coeffcients of the model.
summary(lm_ab)

# Ex 1.4
modelXY <- function(data){
  return(summary(lm(data$Y ~ data$X)))
}

modelXY(GenXY(a = 3,b = 2)) #p value is 0.18 so we do not reject H_0

modelXY(GenXY(a = 3,b = 2, sd = 100)) #p value is 0.0.948 so we do not reject H_0

modelXY(GenXY(a = 3,b = 2, n = 1000)) #p value is very small so we reject H_0

#Ex 1.5 Fit now a regression model without intercept. Compute AIC, BIC for the models with and without intercept. Perform the F-test. Comment the results.
XY <- GenXY(a = 3,b = 2)
intercept <- lm(Y ~ X, data = XY)
no_intercept <- lm(Y ~ X - 1, data = XY)

AIC(intercept, no_intercept)
BIC(intercept, no_intercept)

anova(intercept, no_intercept, test = "F")

## P value is greater than 0.05 so we don't reject H_0 but WHAT IS THE NULL HYPOTHESIS

#Ex 2
#Ex2.1 Generate artificial data from the model
X <- rnorm(50, mean = 0, sd = 1)
Y <- sapply(X, function(x) rnorm(1, mean= x^2 - x + 1, sd = 2))
plot(X, Y)
curve(x^2 -x + 1, add = TRUE, col = "red")

#Ex2.2 Fit a simple linear regression model E(Y jX = x) = beta_0 + beta_1x to the data generated in 2.1. Plot the fitted line on top of the scatter plot as usual.
poly_lm_simple <- lm(Y ~ X)
abline(poly_lm_simple, col = "blue")

#Ex 2.3 Plot the predictor variable vs the residuals and the Q-Q plot of the residuals vs the normal quantiles (qqnorm and qqline functions), comment the plots.
plot(X, residuals(poly_lm_simple))
#Residuals seem random so not correlation with x. This means our distribution is linear
qqnorm(residuals(poly_lm_simple))
qqline(residuals(poly_lm_simple), col = "red")
#There is some deivation at the tails but overall pur samples match the theoretical normalquantiles so it is a normal distribution

#Ex 2.4Fit now the true degree 2 polynomial model (remember that it is still a linear model and we can use the function lm). Plot the result in the same graph with a dierent color, you can also add a legend.
curve(x^2 -x + 1, add = TRUE, col = "red")
abline(poly_lm_simple, col = "blue")
poly_lm_advanced <- lm(Y ~ I(X^2) + X)
abline(poly_lm_advanced, col = "green")
legend("topright", c("True function", "Simple regression", "Advanced regression"), col = c("red", "blue", "green"), lty = 1)

#Ex 2.5 As in 2.3 obtain also for the polynomial regression model the predictor residuals plot and the Q-Q plot vs the normal quantiles.
plot(X, residuals(poly_lm_advanced))

qqnorm(residuals(poly_lm_advanced))
qqline(residuals(poly_lm_advanced), col = "red")

#Ex 2.6 Perform model selection between the simple linear regression of point 2.2 and the polynomial regression in point 2.4. Use the log-likelihood ratio test, the F-test (both with anova). Moreover perform model selection also using AIC and BIC score (AIC, BIC)
AIC(poly_lm_simple, poly_lm_advanced)
BIC(poly_lm_simple, poly_lm_advanced)

anova(poly_lm_simple, poly_lm_advanced, test = "LRT")
anova(poly_lm_simple, poly_lm_advanced, test = "F")
#The advanced model scored lower with both AIC and BIC making it the better choice. Not sure how to interpret anova results

#Ex 2.7 Try to fit now a polynomial of higher degree (e.g. 3,4,5,...). Perform also here model selection. In particular plot the AIC (or BIC) score as a function of the polynomial degree. Plot also the log-likelihood as a function of the polynomial degree. What can you observe? What is the difference between the pure log-likelihood and the AIC and BIC scores?
poly_lm_list <- lapply(1:10, function(i) lm(Y ~ I(X^i)))
poly_lm_list <- lapply(1:10, function(i) lm(Y ~ poly(X, degree = i)))
poly_lm_AIC <- sapply(poly_lm_list, function(x) AIC(x))
poly_lm_LL <- sapply(poly_lm_list, function(x) logLik(x))

plot(1:10, poly_lm_AIC, col = "red", type = "l")
plot(1:10, poly_lm_LL, col = "blue", type = "l")

#The graphs are almost identical in terms of shape. For both graphs the minimum value is at X^3. The odd X^i values are lower then the even X^i values. But this doessnt make any sense since the data was generated with X^2

#Ex 3
#Ex 3.1 As in exercise 1 and 2 generate n = 50 observations from the following model
X <- rnorm(50, mean=0, sd = 0.5) #have to square root variance
Y <- sapply(X, function(x) rnorm(1, mean= exp(-3*x) + 2*x, sd = sqrt(2)))
plot(X, Y)

#Ex 3.2 Fit a simple linear model E(Y jX = x) = beta_0 +beta_1x, and polynomial regression models up to degree 5.
exp_lm_list <- lapply(1:5, function(i) lm(Y ~ I(X^i))) #when X^1 and X are both present, only X^1 is used

#Ex 3.3 Perform model selection of the previous models using AIC and BIC
exp_lm_AIC <- sapply(exp_lm_list, function(x) AIC(x))
exp_lm_BIC <- sapply(exp_lm_list, function(x) BIC(x))
exp_lm_df <- data.frame("Polynomial" = 1:5, "AIC" = exp_lm_AIC, "BIC" = exp_lm_BIC)
min(exp_lm_BIC) # x^2 gives the best model

#Ex 3.4 Check the residuals distribution and the plot the predictor observations vs the residuals. Comment.
plot(1, type = "n", xlim=c(-0.7, 0.7), ylim = c(-4.5, 4.5), xlab = "X", ylab= "Y")
for(i in 1:5){
  points(X, residuals(exp_lm_list[[i]]), col = i+1)
  curve(x^i, add = TRUE, col = i+1)
}
legend("topleft", c("X^1", "X^2", "X^3", "X^4", "X^5"), col = c(2:6), lty = 1, cex = 0.75)
#The higher the power of X, the higher the residuals, so X^1 has the best fit

#Ex 3.5 Now fit the true model E(Y jX = x) = Beta_0 +Beta_1x+exp(Beta_2x) with Gaussian noise.
exp_lm <- nls(Y ~ b0 + b1 * X + exp(b2 * X))
summary(exp_lm)
#nls gives acceptrable results. b1 and b2 are 2 and -3, and nls predicts 2.8 and -2.1 respectively

RSS <- function(parameters, xvals, yvals){
  Y_est <- parameters[1] + parameters[2] * xvals + exp(parameters[3]* xvals)
  square <- (yvals - Y_est)^2
  return(sum(square))
}

optim(par = c(1, 1, -1), fn = RSS, xvals = X, yvals = Y)
#optim of RSS returns results identical to nls

LogLikelihood <- function(parameters, xvals){
  Y_est <- parameters[1] + parameters[2] * xvals + exp(parameters[3]* xvals)
  return(-sum(log(Y_est)))
}

optim(par = c(1, 1, -1), fn = LogLikelihood, xvals = X)
#This gives very bad results. results are better if you remvethe minus sign in the Loglikleihood function, but still not good. I thinj maybe the point of this question is that this method is bad, but I'm not sure

#Ex 4
wines <- read.csv("winequality-red.csv", sep =";")

#Ex 4.1 Fit a linear regression model using all the regressors. Use the function summary, based on the results of the t-test which are the important regressors?
wine_model <- lm(quality ~ ., data = wines)
summary(wine_model)

#The important regressors are alcohol, sulphates, total sulfur dioxide, chlorides and volatile aciditity

#Ex 4.2 Use forward stepwise selection with the AIC score to select the relevant covariates.
wine_model_forward <- lm(quality ~ 1, data=wines)
wine_model_forward
fitAll <- lm(quality ~., data=wines)
wine_model_forward_step <- step(wine_model_forward, scope= formula(fitAll), direction = "forward")

#USED A DIFFERENT METHOD, THIS WORKS BUT ISNT WHAT GHERARDO WANTS, I DONT KNOW HOW TO USE UPDATE
##################
wine_model

lapply(myvars, function(dvar) lm(eval(paste0(dvar,' ~ wt')), data = Boston))
lapply(wines, function(column) update(wine_model_forward, eval(paste0('. ~ ', column)), data = wines))

wine_model_forward <- lm(quality ~ 1, data = wines)

update(wine_model_forward, . ~ wines$chlorides)
lapply(wines, function(x) x/2)

lapply(wines, update(wine_model_forward, . ~ ), data = x))
update(wine_model_forward, wines$chlorides)


#Ex 4.3
st.errors <- summary(wine_model)$coefficients[-1,2]
W <- wine_model$coefficients[-1] / st.errors
ix <- sort(abs(W), decreasing = TRUE, index.return = TRUE)$ix
reg.names <- names(wine_model$coefficients[-1])[ix] ### sorted
sigma2_est <- sum(wine_model$residuals^2) / nrow(wines)

s <- Inf
for (j in 1:length(reg.names)){
  fit.tmp <- lm(paste("quality ~", paste(reg.names[1 : j], collapse = "+") ), data = wines)
  s.tmp <- sum(fit.tmp$residuals^2) + j * (sigma2_est) * log(nrow(wines))
  if (s.tmp < s ){
    J <- j
    s <- s.tmp
  }
}
fit.final <- lm(paste("quality ~", paste(reg.names[1 : J], collapse = "+") ),
                data = wines)
summary(fit.final)


#Ex 4.4

#Ex 5
wines <- read.csv("winequality-red.csv", sep =";")
good <- wines$quality > 5
wines$quality <- "bad"
wines[good, "quality"] <- "good"
wines[,"quality"] <- as.factor(wines[, "quality"])

#Ex 5.1
wine_model_logistic <- glm(quality ~ ., data = wines, family = "binomial")
wine_model_logistic

#Ex5.2
#I tried to find the minimum AIC column in each iteration and then remove that column from the dataframe
wine_quality <- wines #Made this to preserve original dataset
wine_null <- wines #I need one column that doesn't have the quality column otherwise AIC doesnt work
wine_null$quality <- NULL
start_model_wine <- glm(quality ~ 1, family = "binomial", data = wine_quality)

for (i in 1:length(wines)){
  wine_start_AIC <- AIC(start_model_wine)
  AIC_vector <- sapply(wine_null, function(x) AIC(glm(quality ~ x, data = wine_quality, family="binomial"))) #creates list of AICs
  if(min(AIC_vector) < wine_start_AIC){
    start_model_wine <- update(start_model_wine, . ~ . + wines[,match(min(AIC_vector), AIC_vector)]) #updates model with minimum AIC value
    wine_null[,match(min(AIC_vector), AIC_vector)] <- NULL  #removes the minimum column from both dataframes
    wine_quality[,match(min(AIC_vector), AIC_vector)] <- NULL
  }
}
start_model_wine

#hopefully someone can figure this out

####Gherardo
fit <- glm(quality ~ 1, family = "binomial",
           data = wines) ## only the intercept
regressors <- colnames(wines)[-12]
selected <- c()
score <- AIC(fit)
score.best <- score
done <- FALSE
while (!done){
  for (reg in regressors[!(regressors %in% selected)]){
    tmp <- update(fit, formula = paste(". ~ . + ", reg))
    score.tmp <- AIC(tmp)
    if (score.tmp < score.best){
      score.best <- score.tmp
      best <- tmp
      selected.best <- c(selected, reg)
    }
  }
  if (score.best < score){
    fit <- best
    score <- score.best
    selected <- selected.best
  }else{ ### if there is no increase
    done <- TRUE
  }
}


#Ex5.4
preds <- predict(fit)
linkinv <- binomial()$linkinv
linkinv(preds)[1]


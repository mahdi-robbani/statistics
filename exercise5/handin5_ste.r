###Ex 1 Simple linear regression

#Ex 1.1 
sigma = 1
sigma2 = 2
a = 2
b = 0

          # we can write it in two ways
x <- rnorm(50, mean = 0, sd = sigma)   
y <- rnorm(50, mean = a * x + b, sd = sigma2)

          # or
x <- rnorm(50, mean = 0, sd = sigma)     # a = b1, b = b0, rnorm(50) = independent noise
y <- b + a * x + rnorm(50, sd = sigma2)  # y = b + ax + error

plot(x,y)
          # the regression line (curve) consists of the expected values of a variable (Y) 
          # when given the values of independent variable (X). 
          # It is defined as E[Y|X = x]. 
          
curve(b + a * x , add = TRUE, col = "red")  # I guess the true regression function is without the error

# Ex 1.2
model <- lm(y ~ x)
abline(model, col = "blue")  # Plot the ???tted regression line

# Ex 1.3
summary(model)               # Ho : Bo = 0   (intercept = 0)
                             # Ho : B1 = 0   (slope = 0)
                             # if b = 0, p-value > alpha, I don't reject the null
                             # if b = 2, p-value < alpha, I reject the null

# Ex 1.5
model_no.intercept <- lm(y ~ x - 1) # regression model without intercept
summary(model_no.intercept)  

AIC(model, model_no.intercept)
BIC(model, model_no.intercept)      # the model without the intercept is better since it has a lower score

anova(model_no.intercept, model)    # anova(lm.1,lm.2), it performs the F test (if not specified: test = "LRT" etc) to compare lm.1 and lm.2 
                                    # (i.e. it tests whether reduction in the residual sum of squares are statistically significant or not). 
                                    # Note that this makes sense only if lm.1 and lm.2 are nested models.
                                    # if p-value > 0.5, I don't reject the Ho : the two model are the same
                                    # so in this case the simpler model is better
    # in other words a large p-value (larger than alpha) means that adding the intercept did not lead to a 
    # significantly improved fit over the smaller model

    # it looks like we can find the RSS (residual sum of squares) also in this way
deviance(model)
sum(resid(model)^2)


###Ex 2 Poly linear regression

#Ex 2.1
x <- rnorm(50, mean = 0, sd = 1)        # when it's written X ~ N(0,4), in r should we write sd = 2 or 4?
y <- 1 - x + (x^2) + rnorm(50, sd = 2)
plot(x,y)

curve(1 - x + (x^2), add = TRUE, col = "red") # true regression function

#Ex 2.2
model <- lm(y ~ x)
abline(model, col = "blue")
summary(model)

#Ex 2.3
model_res <- resid(model)
plot(x, model_res, xlab = "predictors", ylab = "residuals", main = "Redisuals VS Predictors")
abline(0,0, col = "red")
                    # from the plot we can see that the residuals are homoscedastic, since they 
                    # are randomly dispersed around the horizontal line. Homoscedasticity of
                    # the residuals is one of the conditions we need for linear regression

 # In residual plots, we are looking for the absence of pattern! Any type of pattern exhibited in a residual plot
 # indicates a problem with the model, typically either due to lack of fit or variance heterogeneity.

qqnorm(model_res, ylab = "residuals", xlab = "normal scores", 
       main = " Q-Q plot of the residuals vs the normal quantiles")
qqline(model_res, col = "red")
  # from the QQ plot (Normal Probability Plot of Residuals), we can assess if the error term 
  # is normally distributed (an other condition for the linear regression)

  # Gherardo did not ask for the following plots, but it can be useful to know
par(mfrow = c(2, 2))                          
plot(model)
par(mfrow = c(1, 1))

#Ex 2.4
plot(x,y)         
curve(1 - x + (x^2), add = TRUE, col = "red")
abline(model, col = "blue")                         # 1) Simple model plot

model_poly <- lm(y ~ x + I(x^2))  
points(x, predict(model_poly), col = "orange")      # 2) Poly model plot
summary(model_poly)

model_poly2 <- lm(y ~ poly(x,2))        # it should be the same as the previous one
points(x, predict(model_poly2), col = "dark green")
summary(model_poly2)

#Ex 2.5
plot(x, residuals(model_poly),          # residuals vs predictors plot
     xlab = "predictors", ylab = "residuals", main = "Redisuals VS Predictors")
abline(0,0, col = "red")
                                        # QQ-norm residuals plot
qqnorm(residuals(model_poly), ylab = "residuals", xlab = "normal scores", 
       main = " Q-Q plot of the residuals vs the normal quantiles")
qqline(model_res, col = "red")

#Ex 2.6                                 # LRT, F-test, AIC, BIC
anova(model, model_poly, test = "LRT")
anova(model, model_poly, test = "F")    # p-value is really low, so we reject the Ho : there is
AIC(model, model_poly)                  # not difference between the simple and the complex models
BIC(model, model_poly)                  # (so we should use the complex model)

#Ex 2.7 
model_poly3 <- lm(y ~ x + I(x^2) + I(x^3))
model_poly4 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
model_poly5 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5))
model_poly6 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6))


anova(model_poly, model_poly3, test = "LRT") # p-value > alpha, we don't reject the Ho,   
anova(model_poly, model_poly3, test = "F")   # is not difference in the two model, we should use
                                             # the simpler one


AIC_BIC_LRT_poly <- function (x, y) {        # plot the AIC, BIC and LRT as function of poly degree
  degree <- c()
  BIC <- c()
  AIC <- c()
  LRT <- c()
  for(i in 1:15){
    model_poly_f <- lm(y ~ poly(x,i)) 
    degree <- c(degree, i)
    AIC <- c(AIC, AIC(model_poly_f))
    BIC <- c(BIC, BIC(model_poly_f))
    RSS <- anova(model_poly, model_poly_f, test = "LRT")$RSS[2]
    LRT <- c(LRT, RSS)
  }
  
  output <- list()
  output$degree <- degree
  output$AIC <- AIC
  output$BIC <- BIC
  output$LRT <- LRT
  return (output)
}

AIC_BIC_LRT <- AIC_BIC_LRT_poly(x,y)
plot(AIC_BIC_LRT$degree, AIC_BIC_LRT$AIC, 
     type = "line", col = "red", ylim = c(130, 260),
     main = "AIC, BIC, LRT as function of polynomial degree", 
     ylab = "AIC or BIC scores, RSS for LRT", xlab = "polynomial degree" )
lines(AIC_BIC_LRT$degree, AIC_BIC_LRT$BIC, 
      type = "line", col = "blue", add = TRUE)
lines(AIC_BIC_LRT$degree, AIC_BIC_LRT$LRT, 
      type = "line", col = "dark green", add = TRUE)
legend("bottomleft",legend = c("AIC", "BIC", "LRT"), 
       col = c("blue", "red", "dark green"), lty = 1)

###Ex 3 Non linear regression

#Ex 3.1
x <- rnorm(50, mean = 0, sd = 0.25)
y <- exp(-3*x) + 2*x + rnorm(50, sd = 2)

#Ex 3.2
simple_model <- lm(y ~ x)
model_poly2 <- lm(y ~ x + poly(x,2))
model_poly3 <- lm(y ~ x + poly(x,3))
model_poly4 <- lm(y ~ x + poly(x,4))
model_poly5 <- lm(y ~ x + poly(x,5)) 

#Ex 3.3
AIC(simple_model, model_poly2, model_poly3, model_poly4, model_poly5)
BIC(simple_model, model_poly2, model_poly3, model_poly4, model_poly5) 
     # they have similar score, between these I would use the poly model with 2 degree

#Ex 3.4
model_res <- resid(model_poly2)
qqnorm(model_res, ylab = "residuals", xlab = "normal scores", 
       main = " Q-Q plot of the residuals vs the normal quantiles")
qqline(model_res, col = "red")

plot(x, model_res, xlab = "predictors", ylab = "residuals", main = "Redisuals VS Predictors")
abline(0,0, col = "red")
     # comment? not sure of anything..  

#Ex 3.5
     # wtf??


###Ex 4 Wine quality

#Ex 4.1
wines <- read.csv("winequality-red.csv", sep =";")
model <- lm(quality ~ ., data = wines)
summary(model)
      # the most important regressor should be: alcohol, volatile.acidity, sulphates, 
      # chlorides, total.sulfur.dioxide 

#Ex 4.2
?update
update(model) # ?????

#Ex 4.3
              # ?????



###Ex 5 Logistic regression

              # W T F


###Ex 6

              
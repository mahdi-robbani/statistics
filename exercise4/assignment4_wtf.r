##Ex 1
angles <- read.table("angles.txt")
angles <- as.vector(angles$x)

  # k estimation by MLE numerical optimization from previous exercise
dsink <- function(x, k = 1, lg = FALSE) {
  sinintegral <- integrate(function(x) sin(x) ^ k, lower = 0, upper = pi)$value
  if (lg == FALSE) {
    return(sin(x) ^ k / sinintegral)
  }
  else{
    return(log(sin(x) ^ k / sinintegral))
  }
}

minusll <- function(k, xvals) {
  return(-sum(dsink(xvals, k, lg = TRUE)))
}

k_est <- optimize(f = minusll, xvals = angles, interval = c(0, 100))$minimum

#Ex 1.1 Build 99% CI for k

    #first find SE_est(k) -> bootstrap
vect_k_est_bt <- replicate(1000, expr = { 
  angles_bt <- sample(angles, size = length(angles), replace = TRUE)
  optimize(f = minusll, xvals = angles_bt, interval = c(0, 100))$minimum
})

se_est <- sd(vect_k_est_bt)

    #calculate 99% confidence interval
alpha <- 0.01
z <- qnorm(1 - alpha/2)
a <- k_est - se_est * z
b <- k_est + se_est * z
paste("99% Confidence Interval for k: ( a =", a, ", b= ", b, ")")

#Ex 1.2 Test if k > 10 at a con???dence level ?? = 0.05

alpha <- 0.05                 # Ho = Theta (less or equal) 10
z <- qnorm(1 - alpha)         # H1 = Theta > 10
w <- (k_est - 10) / se_est    
w > z                         # W > Z  TRUE, so we reject Ho and we thinK Theta > 10
                              # with 95% confidence (we can say it because it's
                              # 1 side test, otherwise we will just say that the Ho is
                              # wrong)
qnorm(w)

##Ex.2
spikes <- read.table("neuronspikes.txt")
spikes <- spikes$V1

#2.1
curve(dexp(x, rate = 1), col = "red")
curve(dgamma(x, shape = 1), col = "blue", add = TRUE)
              # skipped 2.2
  
##Ex.3
  # mll exponential
mll_exp <- function(par, xvals){
  return(-sum(dexp(xvals, rate = par, log = TRUE)))
}
exp_par <- optimize(f = mll_exp, interval = c(0,100), xvals = spikes)$minimum


  # mll gamma
mll_gamma <- function(par, xvals){
  return(-sum(dgamma(xvals, shape = par[1], rate = par[2], log = TRUE)))
}
gamma_par <- optim(f = mll_gamma, par = c(1, 1), xvals = spikes)$par

  # mll inverse Gaussian
dinvnorm <- function(x, mu, lambda, lg = FALSE){ 
  if(lg == TRUE){ 
    return(log(sqrt(lambda/(2*pi*(x^3)))*exp(-lambda*((x-mu)^2)/(2*(mu^2)*x)))) 
    } 
  else{ 
    return(sqrt(lambda/(2*pi*(x^3)))*exp(-lambda*((x-mu)^2)/(2*(mu^2)*x))) 
    } 
}
mll_invnorm <- function(par, xvals){ 
  return(-sum(dinvnorm(xvals, mu = par[1], lambda = par[2], lg = TRUE))) 
  }
invnorm_par <- optim(par = c(1, 1), fn = mll_invnorm, xvals = spikes)$par

  # mll log normal
mll_lnorm <- function(par, xvals){
  return(-sum(dlnorm(xvals, meanlog = par[1], sdlog = par[2], log = TRUE)))
}
lnorm_par <- optim(par = c(1, 1), fn = mll_lnorm, xvals = spikes)$par
        # skipped to 4

##Ex.4

#Ex 4.1
ramp_spikes <- read.csv("cell_types.csv", na.strings = "")
ramp_spikes <- ramp_spikes$ef__peak_t_ramp
ramp_spikes <- ramp_spikes[!is.na(ramp_spikes)]

  # mll log normal distribution and MLE
mll_lnormal <- function(par, xvals){
  return(-sum(dlnorm(xvals, meanlog = par[1], sdlog = par[2], log = TRUE)))
}
mu_lnorm_est <-  optim(f = mll_lnormal, par = c(1, 1), xvals = ramp_spikes)$par[1]

  # calculate se of mu estimator (logmu)
vect_lnorm_par_est_bt <- replicate(5, expr = { 
  ramp_spikes_bt <- sample(ramp_spikes, size = length(ramp_spikes), replace = TRUE)
  optim(par = c(1, 1), f = mll_lnormal, xvals = ramp_spikes_bt)$par
})

se_mu_lnorm_est <- sd(vect_lnorm_par_est_bt[1,])

#Ex 4.2 obtain 95% confidence interval for mu
alpha <- 0.05
quantile(vect_lnorm_par_est_bt[1], probs = c(alpha/2, 1 - alpha / 2))

##Ex 5.

#Ex 5.1
ramp_spikes <- read.csv("cell_types.csv")
colnames(ramp_spikes)
ramp_spikes_mouse <- ramp_spikes$ef__peak_t_ramp[ramp_spikes$donor__species == "Mus musculus"]
ramp_spikes_mouse <- ramp_spikes_mouse[!is.na(ramp_spikes_mouse)]
ramp_spikes_human <- ramp_spikes$ef__peak_t_ramp[ramp_spikes$donor__species == "Homo Sapiens"]
ramp_spikes_human <- ramp_spikes_human[!is.na(ramp_spikes_human)]

ramp_spikes_human_log <- log(ramp_spikes_human)
ramp_spikes_mouse_log <- log(ramp_spikes_mouse)

t.test(ramp_spikes_human_log, ramp_spikes_mouse_log) # p-value is < alpha so we reject Ho,
                                                     # Ho: mean_human = mean_mouse
                                                     # H1: mean_
                                                     # so the 2 means are significally different with a 
                                                     # confidence of 95%

  # basic t-test (var.equal = T) works only with samples that have homogeneous variance but by 
  # default the var.equal is set to FALSE and so we can use it with sample that have different variance
  
#Ex 5.2
alpha <- 0.05
Z005 <- qnorm(1 - alpha/2)

Zobserved <- (mean(ramp_spikes_human_log) - mean(ramp_spikes_mouse_log)) / (
  sqrt(
    (var(ramp_spikes_human_log)^2 / length(ramp_spikes_human_log)) + 
      (var(ramp_spikes_mouse_log)^2/length(ramp_spikes_mouse_log)
      )
  ))
Zobserved <= Z005        # Ho: mean_human = mean_mouse
                         # since Z (or W) rejection of Ho
                         # two sample wald test? looks like we did not cover 
                         # it but this should be correct

##Ex 6

#Ex 6.1
x = rnorm(20, mean = 2, sd = 4)
y = rnorm(40, mean = 2.5, sd = 4) 

#Ex 6.2 
Txy <- (sqrt((length(x)+length(y))/length(x)*length(y))*(mean(x)-mean(y)))/sd(x+y)
p_value <- 2*pt(q = -abs(Txy), df = length(x)+length(y)-2)
paste(p_value)

t.test(x, y, var.equal = TRUE)    # there may be a problem since we expected the p-value 
                                  # from t.test and from analytical formula to be the same
                                  # but it looks like they are not

#Ex 6.3                           # we miss Wald test here
vect_norm_meandiff_bt <- replicate(10000, expr = {
  x_bt <- sample(x, size = length(x), replace = TRUE)
  y_bt <- sample(y, size = length(y), replace = TRUE)
  mean(x_bt)-mean(y_bt)
})
se_meandiff <- sd(vect_norm_meandiff_bt)

#Ex 6.4

#Ex 6.5
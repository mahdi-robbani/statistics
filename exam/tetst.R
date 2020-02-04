S <- rnorm(1000, mean = -2, sd = 5)

mu_est <- mean(S)


sd_est <- sqrt( var(S)*(length(S) - 1)/ length(S) )

paste("mean:", mu_est, " sd:", sd_est)

#numerical
mllk <- function(pars, data){
  if (pars[2] <= 0){
    return( Inf)  ## sigma can not be negative
  }
  sum( -dnorm(x = data, mean = pars[1], sd = pars[2], log = T))
}

optim(par = c(0, 1), fn = mllk, data = S)

#derivative based

gmllk <- function(pars, data){
  dm <- -sum((data - pars[1])) / (pars[2] ^ 2)
  ds <- -sum((data - pars[1]) ^ 2) / (pars[2] ^ 3) + length(data) / (pars[2])
  return(c(dm, ds))
}

optim(par = c(3,8), fn = mllk, gr = gmllk, data = S, method = "BFGS")

optim(par = c(3,8), fn = mllk, data = S, method = "BFGS")

#gradient descent
#### plotting the contour plot
#mu_grid <- seq(-9.9, 9.9, 0.2)
#sigma_grid <- seq(0, 9.8, 1)
mu_grid <- seq(from = -10, to = 10, length.out = 100)
sigma_grid <- seq(from = 1, to = 10, length.out = 100  )
grid <- expand.grid(mu_grid, sigma_grid)

z <- apply(grid, MARGIN = 1, mllk, data = S)
Z <- matrix(data = z, byrow = F, nrow = 100)
contour(x = mu_grid, y = sigma_grid, z = Z, xlab = "mu", ylab = "sigma",
        levels = quantile(z, probs = seq(from = 0, to =1, by = 0.1)),
        main = paste0("Contour plot of  minus log likelihood"))
par0 <- c(3, 8) ##initial parameters guess
k <- 0.01 ## learning rate of gradient descent
M <- 40 ## max num of iterations
for (i in 1:M){
  parOld <- par0
  par0 <- par0 - k*gmllk(par0, S)  ##gradient descent update
  points(par0[1], par0[2], col = "red", lwd = 1)
  arrows(x0 = parOld[1], y0 = parOld[2] , x1 = par0[1], y1 = par0[2],
         col = "red",length = 0.01, lwd = 1 )
}

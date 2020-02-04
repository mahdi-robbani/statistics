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
  JK_pairs[[i]] <- simulate_JKpair(t=10, alpha = 1) 
}
###
tval <- 1
aval <- 0.1
JK_pairs <- list()
for (i in 1:1000){
  JK_pairs[[i]] <- simulate_JKpair(start = sample(c("A", "C", "T", "G"), size = 1),
                                   t=tval, alpha = aval) 
}



#Minus log likelihood function
distJKMLL <- function(alpha, xvals, t = 10){
  probs <- sapply(xvals, function(x) JKprob(X = x[1], Y = x[2], alpha = alpha, t = t))
  return(-sum(log(probs)))
}

#Solving numerically
optimize(f = distJKMLL, xvals = JK_pairs, interval = c(0, 0.01), t = 10)


############
####weekly exercise solution
JK_table <- data.frame(x = sapply(JK_pairs,function(x) x[1]), y= sapply(JK_pairs,function(x) x[2]))

f_jk <- function(x, y, t, a){
  if (x == y){
    0.25 + 0.75 * exp(-4*a*t)
  }else{
    0.25 - 0.25 * exp(-4*a*t)
  }
}


mll_jk <- function(a, data){
  -sum(sapply(1:nrow(data), function(i){
    log(f_jk(x = data$x[i], y = data$y[i], t = tval, a = a))
  }))
}


a_est <- optimize(f = mll_jk, interval = c(0,10), data = JK_table, t = tval)$minimum

a_est


##analytical
n1 <- sum(apply(JK_table, MARGIN = 1,
                function(r) r[1] == r[2]))
n2 <- nrow(JK_table) - n1
a_mle <- -log( (3*n1 - n2)/ (3*(n1+n2))) / (4*tval)
a_mle

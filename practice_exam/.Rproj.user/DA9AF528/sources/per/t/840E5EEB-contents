

dsin <- function(x, k = 1){
  const <- integrate(function(t){ sin(t)^k }, lower = 0, upper = pi)$value
  if (x < 0 || x > pi){
    return(0)
  }
  return( (sin(x) ^ k )/ const )
}

psin <- Vectorize(function(q, k = 1){
  if (q > pi){
    return(1)
  }
  if ( q < 0){
    return(0)
  }
  const <- integrate(function(t){ sin(t)^k }, lower = 0, upper = pi)$value
  integrate(function(t){ sin(t)^k }, lower = 0, upper = q)$value / const
}, vectorize.args = "q")

### to be vectorized 
qsin <- function(p, k = 1){
  uniroot(function(x){
    psin(x, k) - p   
  }, interval = c(0, pi), extendInt = "upX")$root
}

rsin <- function(n, k = 1) {
  return(sapply(runif(n), qsin, k = k))
}
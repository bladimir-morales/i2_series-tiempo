## Durbin-Levinson:
DurbinLevinson = function(serie, ar = NULL, ma = NULL){
  X = serie
  nu = c()# vector de nu
  X.hat = c()## vector de predicciones a un paso
  N = length(X)## largo del vector
  Phi = matrix(0, ncol = N, nrow = N)## Matriz de coeficientes phi[i,j]
  rho = ARMAacf(ar = ar, ma = ma, lag.max = N)[2:(N+1)]
  Phi[1,1] = rho[1]
  nu[1] = (1 - Phi[1,1]^2)
  X.hat[1] = 0
  X.hat[2] = Phi[1,1]*X[1]
  for(n in 2:(N-1)){
    Phi[n,n] = (rho[n] - sum(Phi[n - 1,1:(n - 1)]*rho[n - 1:(n - 1)]))/nu[n - 1]
    Phi[n,1:(n-1)] = Phi[n - 1,1:(n - 1)] - Phi[n,n]*Phi[n - 1,(n - 1):1]
    nu[n] = nu[n-1]*(1-Phi[n,n]^2)
    X.hat[n+1] = sum(X[n:1]*Phi[n,1:n])
  }
  nu = ((N-1)*var(serie)/N)*c(1,nu)
  list(fitted = X.hat, nu = nu)
}
## Durbil-Levinson log likelihood
dl.loglik = function(x, serie, p = 1, q = 1){
  ar = NULL
  ma = NULL 
  if(p > 0){
    ar  = x[1:p]
  }
  if(q > 0){
    ma  = x[(p+1):(p+q)]
  }
  fit   = DurbinLevinson(serie = serie, ar = ar, ma = ma)
  e     = serie - fit$fitted
  nu    = fit$nu 
  aux = 0.5*(sum(log(nu)) +  sum(e^2/nu))
  if(is.nan(aux)){aux = Inf}
  aux
}

# Bartlett
Bartlett <- function(ar = phi, lag.max = 10){
  N = 100
  acf_teo = ARMAacf(ar = phi, lag.max = lag.max + N)
  barlett = rep(0,lag.max)
  for (h in 1:lag.max){
    for (j in h:N){
      barlett[h] = barlett[h] + (acf_teo[(j + h) + 1] + 
                                   acf_teo[(j - h) + 1] - 
                                   2*acf_teo[j + 1]*acf_teo[h + 1])^2
    }
  }
  return (barlett)
}

# Whittle
w.loglik <- function(x, serie, p = 1, q = 1){
  X <- serie
  if(p > 0){
    ar <- x[1:p]
  }
  if(q > 0){
    ma <- x[(p+1):(p+q)]
  }
  if(p == 0){ar = numeric()}
  if(q == 0){ma = numeric()}	
  sigma <- x[p+q+1]
  n <- length(X)
  aux <- LSTS::periodogram(X, plot = F)
  lambda <- aux$lambda
  I <- aux$periodogram
  f <- LSTS::spectral.density(ar = ar, ma = ma, sd = sigma, lambda = aux$lambda)
  aux <- ( sum(log(f))+ sum(I/f))/(2*n)
  aux
}


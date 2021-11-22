MakeDesignMatrix = function(n, p){
  rho = 0.5
  sigma = matrix(rho, p, p)
  diag(sigma) = diag(sigma) + (1-rho)
  L = chol(sigma)
  z = matrix(rnorm(n * p), ncol = p)
  return(scale(z %*% L))
}

MakeResponses = function(X, SNR, k0 = 10){
  n = nrow(X); p = ncol(X)
  epsilon = matrix(rnorm(n), ncol = 1)
  beta0 = matrix(c(rep(1, k0), rep(0, p-k0)), ncol = 1)
  Fmean = X %*% beta0
  var_nf = var(Fmean)*(n-1)/n; var_ne = var(epsilon)*(n-1)/n
  sigma0 = sqrt(var_nf/(SNR * var_ne))[1,1]
  Y = Fmean + sigma0 * epsilon
  return(list(Y = Y, beta0 = beta0))
}

GenerateReplicate = function(n, p, SNR, k0 = 10){
  X = MakeDesignMatrix(n, p)
  response_list = MakeResponses(X, SNR, k0 = k0)
  Y = response_list$Y
  beta0 = response_list$beta0
  return(list(X = X, Y = Y, beta0 = beta0))
}

test = F
if(test){
  set.seed(2)
  n = 10000; p =  200; SNR = 2
  replicate_list = GenerateReplicate(n, p, SNR)
  X = replicate_list$X
  Y = replicate_list$Y
  beta0 = replicate_list$beta0
  plot(X[,1], X[,2])
  plot(Y, X[,1])
  #print(Y)
  hist(Y)
}
  

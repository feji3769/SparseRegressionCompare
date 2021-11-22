horseshoe = function(y,X, tau = 1, Sigma2 = 1,
                     burn = 1000, nmc = 5000, thin = 1, alpha = 0.05, verbose = F)
{

  N=burn+nmc; effsamp=(N-burn)/thin
  n=nrow(X); p=ncol(X)
  
  ## parameters ##
  Beta=rep(0,p); lambda=rep(1,p);
  sigma_sq = Sigma2;
  
  ## output ##
  betaout=matrix(0,p,effsamp)
  tauout=rep(0,effsamp)
  sigmaSqout=rep(0,effsamp)
  
  ## which algo to use ##
  if(p>n)
    algo=1
  else
    algo=2
  ## matrices ##
  I_n=diag(n);l0=rep(0,p)
  l1=rep(1,n); l2=rep(1,p)
  if(algo==2)
  {
    Q_star=t(X)%*%X
  }
  ## start Gibb's sampling ##
  for(i in 1:N)
  {
    ## update beta ##
    if(algo==1)
    {
      lambda_star=tau*lambda
      U=as.numeric(lambda_star^2)*t(X)
      ## step 1 ##
      u=stats::rnorm(l2,l0,lambda_star)
      v=X%*%u + stats::rnorm(n)
      ## step 2 ##
      v_star=solve((X%*%U+I_n),((y/sqrt(sigma_sq))-v))
      Beta=sqrt(sigma_sq)*(u+U%*%v_star)
    }
    else if(algo==2)
    {
      lambda_star=tau*lambda
      # Q=(1/sigma_sq)*(t(X)%*%X+diag(1/lambda_star^2,p,p))#
      # b=(t(y)*X)/sigma_sq #
      L=chol((1/sigma_sq)*(Q_star+diag(1/as.numeric(lambda_star^2),p,p)))
      v=forwardsolve(t(L),t(t(y)%*%X)/sigma_sq) # changed from solve.
      mu=backsolve(L,v) # changed from solve.
      u=backsolve(L,stats::rnorm(p)) # changed from solve.
      Beta=mu+u
    }
    
    
    ## update lambda_j's in a block using slice sampling ##
    eta = 1/(lambda^2)
    upsi = stats::runif(p,0,1/(1+eta))
    tempps = Beta^2/(2*sigma_sq*tau^2)
    ub = (1-upsi)/upsi
    # now sample eta from exp(tempv) truncated between 0 & upsi/(1-upsi)
    Fub = 1 - exp(-tempps*ub) # exp cdf at ub
    Fub[Fub < (1e-4)] = 1e-4;  # for numerical stability
    up = stats::runif(p,0,Fub)
    eta = -log(1-up)/tempps
    lambda = 1/sqrt(eta);
    
    ## update tau ##
    tempt = sum((Beta/lambda)^2)/(2*sigma_sq)
    et = 1/tau^2
    utau = stats::runif(1,0,1/(1+et))
    ubt = (1-utau)/utau
    Fubt = stats::pgamma(ubt,(p+1)/2,scale=1/tempt)
    Fubt = max(Fubt,1e-8) # for numerical stability
    ut = stats::runif(1,0,Fubt)
    et = stats::qgamma(ut,(p+1)/2,scale=1/tempt)
    tau = 1/sqrt(et)

    ## update sigma_sq ##

    if(algo==1)
    {
      E_1=max(t(y-X%*%Beta)%*%(y-X%*%Beta),(1e-10))
      E_2=max(sum(Beta^2/((tau*lambda))^2),(1e-10))
    }
    
    else
    {
      E_1=max(t(y-X%*%Beta)%*%(y-X%*%Beta),1e-8)
      E_2=max(sum(Beta^2/((tau*lambda))^2),1e-8)
    }
    sigma_sq=1/stats::rgamma(1,(n+p)/2,scale=2/(E_1+E_2))

    if ((i%%1000 == 0) & (verbose))
    {
      print(i)
    }
    
    if(i > burn && i%%thin== 0)
    {
      betaout[,(i-burn)/thin] = Beta
      tauout[(i-burn)/thin]=tau
      sigmaSqout[(i-burn)/thin]=sigma_sq
    }
  }
  pMean=apply(betaout,1,mean)
  result=list("BetaHat"=pMean, "BetaSamples"=betaout)
  return(result)
}

HorseshoeMCMC = function(X,Y,nmc = 1000, burn = 500){
  startTime = proc.time()
  InferenceResultList = horseshoe(Y,X, nmc = nmc,burn = burn)
  totalTime = unname((proc.time() - startTime)[3])
  InferenceResultList$CPUTime = totalTime/(nmc + burn)
  return(InferenceResultList)
}

PostHorseshoeMCMC = function(InferenceResultList){
  beta_hats = abs(InferenceResultList$BetaHat)
  kmeans = kmeans(beta_hats, 2, algorithm = c("Lloyd"))
  signal_id = c(1,2)[(kmeans$centers == max(kmeans$centers))]
  return(1 * (kmeans$cluster == signal_id))
}

test_Horseshoe = F
if(test_Horseshoe){
n = 100; p = 15; SNR = 4
#replicate_list = GenerateReplicate(n, p, SNR)
X = replicate_list$X
Y = replicate_list$Y
beta0 = replicate_list$beta0
InferenceResultList = horseshoe(Y, X, nmc = 5000,burn = 500)
plot(abs(InferenceResultList$BetaHat))
}

coeffDetermination = function(X, gmma, y, y_mag){
  S = (gmma == 1); p = ncol(X)
  if(sum(S) > 0){
    X_g = X[,gmma == 1]
    M = t(X_g) %*% X_g
    if (1/kappa(M) > 1e-4){
      phi_g = X_g %*% solve(M) %*% t(X_g)
    } else {
      phi_g = X_g %*% MASS::ginv(M) %*% t(X_g)
    }
    return((t(y) %*% phi_g %*% y) / y_mag)
  }else {
    return(0)
  }
} 
logPrior = function(gmma, kappa, S0){
  S = (gmma > 0);p = length(gmma)
  if(sum(S) > S0){
    return(-100000)
  } else{
    return( -1 * kappa * sum(S) * log(p))
  }
}
logPriorRatio = function(gmma_p, gmma, kappa, S0){
  return(logPrior(gmma_p, kappa, S0)-logPrior(gmma, kappa, S0))
}
logLikelihoodRatio = function(X, g, gmma_p, gmma, y, y_mag){
  n = nrow(X)
  # terms for likelihood for gmma prime.
  pg_p = sum(gmma_p)
  Rp_p = coeffDetermination(X, gmma_p, y, y_mag)
  log_lik_p = -0.5 * pg_p * log(1+g) - n/2 * log(1 + g*(1-Rp_p))
  # terms for likelihood for previous gmma.
  pg = sum(gmma)
  Rp = coeffDetermination(X, gmma, y, y_mag)
  log_lik = -0.5 * pg * log(1+g) - n/2 * log(1 + g*(1-Rp))
  return(log_lik_p - log_lik)
}
logS = function(gmma){
  p = length(gmma)
  if((sum(gmma) != 0) & (sum(gmma) != p)){
    return(log( 0.5/p + 0.5 / (sum(gmma) * sum(1-gmma)) ))
  } else {
    return(-1*log(p))
  }
}
logSratio = function(gmma_p, gmma){
  log_S_gmma_p = logS(gmma_p)
  log_S_gmma = logS(gmma)
  return(log_S_gmma_p - log_S_gmma)
}
proposeN1 = function(gmma){
  i = sample(1:length(gmma),1)
  gmma[i] = 1-gmma[i]
  return(gmma)
}

proposeN2 = function(gmma){
  p = length(gmma); indx = 1:p
  S = indx[gmma > 0]; Sc = indx[gmma == 0];
  nS = sum(gmma); nSc = sum(1-gmma)
  if((nS > 0) & (nSc > 0)){
    iS = sample(S, 1); iSc = sample(Sc, 1)
    gmma[iS] = 1-gmma[iS]; gmma[iSc] = 1-gmma[iSc]
  } else {
    gmma = proposeN1(gmma)
  }
  return(gmma)
}

logAcceptProb = function(X, g, gmma_p, gmma, y, y_mag, kappa = 1.5){
  S0 = nrow(X)
  log_like_rat = logLikelihoodRatio(X, g, gmma_p, gmma, y, y_mag)[1,1]
  log_p_rat = logPriorRatio(gmma_p, gmma, kappa, S0)
  log_S_rat = logSratio(gmma_p, gmma)
  return(min(1, exp(log_p_rat + log_like_rat + log_S_rat)))
}

FlippingMCMC = function(X,Y, nmc = 1000, burn = 500, verbose = F){
  p = ncol(X); g = p ** 2
  y_mag = sum(Y ** 2)
  gmma_samples = matrix(NA, nrow = burn + nmc, ncol = p)
  gmma_samples[1, ] = sample(c(0,1), p, replace = T, prob = c(.8,.2))
  z = runif(burn + nmc)
  u = runif(burn + nmc)
  startTime = proc.time()
  for(i in 2:(burn + nmc)){
    if((i %% 100 == 0) & verbose){
      print(i)
    }
    if(z[i] <= 0.5){
      proposal = proposeN1(gmma_samples[i-1,])
    } else{
      proposal = proposeN2(gmma_samples[i-1,])
    }
    log_R_g_gp = logAcceptProb(X, g, proposal, gmma_samples[i-1,], Y, y_mag)
    if(log_R_g_gp >= u[i]){
      gmma_samples[i,] = proposal
    } else {
      gmma_samples[i,] = gmma_samples[i-1, ]
    }
  }
  totalTime = unname((proc.time() - startTime)[3])
  InferenceResultList = list(
              burn = gmma_samples[1:burn, ], 
              samples = gmma_samples[(burn+1):(burn + nmc),], 
              CPUTime = totalTime/(burn + nmc)
              )
  return(InferenceResultList)
}


PostFlippingMCMC = function(InferenceResultList){
  nmc = nrow(InferenceResultList$samples)
  incl_prob = apply(InferenceResultList$samples, 2, sum)/nmc
  return(1 * (incl_prob > 0.5))
}

test_flipping = F
if(test_flipping){
  set.seed(1)
  n = 100; p = 200; SNR = 4
  replicate_list = GenerateReplicate(n, p, SNR)
  X = replicate_list$X
  Y = replicate_list$Y
  Y = Y
  beta0 = replicate_list$beta0
  g = p ** 2
  num_burn = 1000
  num_iter = 2000
  InferenceResultList = FlippingMCMC(X,Y, nmc = num_iter, burn = num_burn, verbose = T)
  PostFlippingMCMC(InferenceResultList)
}


SpikeSlabVB = function(X, Y, verbose = F){
  startTime = proc.time()
  res = varbvs::varbvs(X, Z = NULL, y = Y, family = "gaussian", 
                 verbose = verbose)
  totalTime = unname((proc.time() - startTime)[3])
  InferenceResultList = list(
    CPUTime = totalTime, 
    pip = res$pip)
  return(InferenceResultList)
}
PostSpikeSlabVB = function(InferenceResultList){
  return(unname(1 * (InferenceResultList$pip > 0.5)))
}

test_varbvs = T
if(test_varbvs){
  set.seed(1)
  n = 100; p = 200; SNR = 4
  replicate_list = GenerateReplicate(n, p, SNR)
  X = replicate_list$X
  Y = replicate_list$Y
  Y = Y
  beta0 = replicate_list$beta0
  g = p ** 2
  num_burn = 1000
  num_iter = 2000
  InferenceResultList = SpikeSlabVB(X,Y,  verbose = F)
  res = PostSpikeSlabVB(InferenceResultList)
}

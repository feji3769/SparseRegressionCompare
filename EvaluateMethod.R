GetZeroOneError = function(beta0_hat, beta0){
  return(sum(!(beta0_hat == beta0)))
}

GetSupportSize = function(beta0_hat){
  return(sum(beta0_hat == 1))
}

RunExperiment = function(n, p, SNR, InferenceMethod,
                         PostProcessMethod,k0 = 10){
  replicate_list = GenerateReplicate(n, p, SNR)
  X = replicate_list$X
  Y = replicate_list$Y
  beta0 = replicate_list$beta0
  InferenceResultList = InferenceMethod(X,Y)
  CPUTime = InferenceResultList$CPUTime
  beta0_hat = PostProcessMethod(InferenceResultList)
  zeroOneError = GetZeroOneError(beta0_hat, beta0)
  supportSize = GetSupportSize(beta0_hat)
  experiment_result = list(
    zeroOneError = zeroOneError, 
    supportSize = supportSize,
    CPUTime = CPUTime
  )
  return(experiment_result)
}


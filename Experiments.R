#!/usr/bin/env Rscript
library(argparse)
p = ArgumentParser(description = "Run comparison of sparse regression models.")
p$add_argument("--r", type = "integer", help = "Replicated ID.", default = 0)
p$add_argument("--n", type = "integer", help = "Number of observations.", default = 100)
p$add_argument( "--p", type = "integer", help = "Dimension of each input.", default = 100)
p$add_argument("--SNR", type = "double", help = "Signal to noise ratio.", default = 4)
p$add_argument("--method", type = "character", 
               help = "Inference technique: MCMCHorseshoe, flippingSpikeSlab, VBSpikeSlab.",
               default = "flippingSpikeSlab")



source("EvaluateMethod.R")
source("MakeData.R")
source("Methods.R")
argv = p$parse_args()

InferenceMethod = switch(
  argv$method, 
  "MCMCHorseshoe"=HorseshoeMCMC,
  "flippingSpikeSlab" = FlippingMCMC, 
  "VBSpikeSlab" = SpikeSlabVB
)
PostProcessMethod = switch(
  argv$method, 
  "MCMCHorseshoe"=PostHorseshoeMCMC,
  "flippingSpikeSlab"=PostFlippingMCMC,
  "VBSpikeSlab" = PostSpikeSlabVB
)

result = RunExperiment(argv$n, argv$p, argv$SNR, InferenceMethod, PostProcessMethod,k0 = 10)
data = data.frame(
  n = c(argv$n), 
  p = c(argv$p), 
  SNR = c(argv$SNR), 
  method = c(argv$method), 
  zeroOneError = result$zeroOneError, 
  supportSize = result$supportSize, 
  CPUTime = unname(result$CPUTime[1]), 
  ReplicateId = argv$r)
append_bool = file.exists("experiment_result.csv")
write.table(x = data, file = "experiment_result.csv", 
          append = append_bool, sep = ",",
          row.names = F, col.names = !append_bool)





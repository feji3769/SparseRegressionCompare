#!/usr/bin/env Rscript
library(argparse)
parser = ArgumentParser(description = "Run comparison of sparse regression models.")
parser$add_argument("-r", type = "integer", help = "Replicated ID.", default = 0)
parser$add_argument("-n", type = "integer", help = "Number of observations.", default = 100)
parser$add_argument("-p", type = "integer", help = "Dimension of each input.", default = 100)
parser$add_argument("-SNR", type = "double", help = "Signal to noise ratio.", default = 4)
parser$add_argument("-method", type = "character", 
               help = "Inference technique: MCMCHorseshoe, flippingSpikeSlab, VBSpikeSlab.",
               default = "flippingSpikeSlab")

argv = parser$parse_args()

source("EvaluateMethod.R")
source("MakeData.R")
source("Methods.R")


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
print(sprintf("Replicate : %d, Method : %s, n : %d, p : %d, SNR : %f", 
argv$r, argv$method, argv$n, argv$p, argv$SNR))
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
append_bool = file.exists("experimentResults/experiment_result.csv")
write.table(x = data, file = "experimentResults/experiment_result.csv", 
          append = append_bool, sep = ",",
          row.names = F, col.names = !append_bool)





#!/usr/bin/env Rscript
## Dani & Robin
## Initialize Sims
#:::::::::::::::::::::::::::::::::::::::::::
args <- commandArgs(trailingOnly = TRUE)
if ( length(args) < 1 ) {
  cat("Need species name command line parameter")
}
fishSpecies <- args[1]
# get species info 
source(paste(fishSpecies, "Species_Info.R", sep ="_"))
source("Model_Dynamics.R")
source("Helper_Functions.R")
# run burnin 
BurnIn_Reps = 50
Collect_Reps = 200
Initialize_List = InitializePops(2)
Burnin_List = RunSenario_Burnin(Initialize_List, NYears_Burnin = BurnIn_Reps + Collect_Reps)
save(Burnin_List, file = paste(fishSpecies, "Initialized_Pops.RData", sep ="_"))

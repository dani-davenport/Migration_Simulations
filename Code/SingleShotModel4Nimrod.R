#!/usr/bin/env Rscript
## Dani & Robin
## Population Dynamics
## Run Simulations
#:::::::::::::::::::::::::::::::::::::::::::::::::
source("Model_Dynamics.R")
source("Helper_Functions.R")
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
if ( length(args) < 3 ) {
  cat("Need 3 command line parameters i.e.  lambda numRep species\n")
}
lambda <- as.numeric(args[1])
numRep <- as.numeric(args[2])
fishSpecies <- args[3]
OUTPUT = paste(fishSpecies, lambda, sep = "_")
#:::::::::::::::::::::::::::::::::::::::::::::::::
## Initialize Pops
# get data generatated 
source(paste(fishSpecies, "Species_Info.R", sep = "_"))
load(paste(fishSpecies, "Initialized_Pops.RData", sep ="_"))
#:::::::::::::::::::::::::::::::::::::::::::::::::
## Run Model
#:::::::::::::::::::::::::::::::::::::::::::::::::
Migrate_Reps = 500
OUTPUT<- RunSenario_Migrate(Burnin_List, lambda = lambda , NYears_Migrate = Migrate_Reps)

save(list = ls(all.names = TRUE), file = paste(fishSpecies, paste0(paste(lambda, numRep, sep = "_"), "_Simulation.RData")))

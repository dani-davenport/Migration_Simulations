#!/sw7/RCC/NimrodG/embedded-1.11.0/bin/nimexec
#PBS -N NimrodFish
#PBS -A UQ-Medicine
#PBS -l walltime=168:00:00
#PBS -l select=1:ncpus=4:ompthreads=1:mem=20GB

## use this one for the real thing
## #PBS -l select=5:ncpus=20:ompthreads=1:mem=120GB

#NIM shebang /bin/bash

# Repeats r #FIX THE RANGE
#NIM parameter r label "Repeats" integer range from 1 to 5 step 1

# Lambda l This gives you the actual lambda values you require (as floats)
#NIM parameter l label "Lambda" float select anyof 0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 0.2

#NIM parameter s text select anyof "Hake" "Sablefish" "Pollock" "Herring"

#Make sure you submit the job from the desired directory (one with the R code and Outputs directory)
cd $PBS_O_WORKDIR

#Ensure that the Outputs directory exists
if [ ! -d Outputs ]; then
  mkdir Outputs
fi

# Ensure that the Pops have been initialised
if [ ! -f ${NIMROD_VAR_s}_Initialized_Pops.RData ]; then
  Rscript Initialize_Pops.R $NIMROD_VAR_s
  exit 0
fi

#We want to skip over any that have previously been completed.
#If you know what the output filepath is, then use that.
#In this test job it just touched a blank file.
if [ -f Outputs/${NIMROD_VAR_s}_${NIMROD_VAR_l}_${NIMROD_VAR_r}_Simulation.RData ]; then
  #the output file already exists so we don't repeat ... exit gracefully)
  echo "Outputs/${NIMROD_VAR_l}_${NIMROD_VAR_r}_${NIMROD_VAR_s} already exists."
  exit 0
fi


module load anaconda
source activate /groups/MolFishLab/seqDataPy3.6

#The R code expects the 3 parameters as follows
#lambda <- as.numeric(args[1]) #numRep <- as.numeric(args[2]) #fishSpecies <- args[3]

Rscript SingleShotModel4Nimrod.R $NIMROD_VAR_l  $NIMROD_VAR_r $NIMROD_VAR_s

touch Outputs/${NIMROD_VAR_s}_${NIMROD_VAR_l}_${NIMROD_VAR_r}_Simulation.RData

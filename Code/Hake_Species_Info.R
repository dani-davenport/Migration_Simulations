#!/usr/bin/env Rscript
## Dani & Robin
## Population Dynamics
## Genetic and Demographic Indepedance
#:::::::::::::::::::::::::::::::::::::::::::::::::
# Simulations
#:::::::::::::::::::::::::::::::::::::::::::::::::
# Required Libraries & Functions

#:::::::::::::::::::::::::::::::::::::::::::::::::
## Define Species Parameters
#:::::::::::::::::::::::::::::::::::::::::::::::::
# Hake Info
alpha = 3  # age at maturity, assumed fixed
MaxAge = 21  ## oldest age for reproduction, after which all die
####### For initilaistion
# Make survival variable (in function) - SurvivalVariable
# Set up these to initilize pops
Survival = rep(0.791,MaxAge)  ## constant survival
Lx = cumprod(Survival)/Survival  ## cumulative survival
# Number of recruits is variable (in function), but set here to initialize
N1 = 5000
Nx = round(N1*Lx)  ## stable age composition
NTot = sum(Nx)  #sum(Nx), total population of 1+ individuals
AdultN = sum(Nx[alpha:MaxAge])
Fa = c(0,	0,	0.000,	0.321,	0.466,	0.49,	0.539,	0.598,	0.681,	0.726,	0.837,	0.89,	0.964,	1.015,	0.958,	0.928,	0.928,	0.928,	0.928,	0.928,	0.928)  ## relative, age-specific fecundity, length needs to match MaxAge
#                   0.067
Bx = Fa*Nx  ## relative number of offspring produced by each age class
FF = Bx/(sum(Bx))   ## F rescaled to sum to 1 = probability an offspring has parent of age x
units = seq(1:MaxAge)
G = sum(units*Fa*Lx)/sum(Fa*Lx)  ## generation length
L = 10 ## number of diallelic loci

# Popultion Dynamics
# Step 1. Scale by bx to produce N1 recruits, this returns the number of Eggs needed
NEggs1 = sum(Bx * Nx) ###  total egg production at stable N
# mortality alpha
## mortality rate required to produce N1 recruits, given nominal bx
M_a= (N1/NEggs1)
M_a
adjust = (N1/NEggs1)/M_a  ## adjust > 1 indicates bx is too low, given M
adjust
# for herring, the value is not greater than one (it is 1), so we do not need to use the Bx adjusted values
bxprime = Bx *adjust
# check
Eggs2 = bxprime*Nx
Teggs2 = sum(Eggs2)  ## eggs required to produce desired number of recruits (N1)
Teggs2*M_a  ## should equal N1

## adjust
MAlpha = M_a # Gompertz alpha
BHBeta = M_a # Becomes 1-B


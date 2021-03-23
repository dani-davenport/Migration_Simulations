#!/usr/bin/env Rscript
## Dani & Robin
## Population Dynamics
## Genetic and Demographic Indepedance 
#:::::::::::::::::::::::::::::::::::::::::::::::::
# Simulations 
#::::::::::::::::::::::::::::::::::::::::::::::::: 
# Required Libraries & Functions


#:::::::::::::::::::::::::::::::::::::::::::::::::
# # Herring Info
alpha = 2  # age at maturity, assumed fixed # Okamoto lists 2 = 0.2, 3 = 0.9, 4+ = 1
MaxAge = 10  ## oldest age for reproduction, after which all die
####### For initilaistion
# Make survival variable (in function) - SurvivalVariable
# Set up these to initilize pops
Survival = rep(0.6,MaxAge)  ## constant survival
Lx = cumprod(Survival)/Survival  ## cumulative survival
# Number of recruits is variable (in function), but set here to initialize
N1 = 10000
Nx = round(N1*Lx)  ## stable age composition
NTot = sum(Nx)  #sum(Nx), total population of 1+ individuals
AdultN = sum(Nx[alpha:MaxAge])
Wx = sapply(list(seq(2:10)), function(a) 4.5*(10^-6)*(27*(1-exp(0.48*a)))^3.127) # weight at age
Fa = c(0,0,2,3,4,5,6,7,8,9)   ## get larger with increasing size for herring, relative, age-specific fecundity, length needs to match MaxAge
Bx = Fa*Nx  ## relative number of offspring produced by each age class
FF = Bx/(sum(Bx))   ## F rescaled to sum to 1 = probability an offspring has parent of age x
units = seq(1:MaxAge)
G = sum(units*Fa*Lx)/sum(Fa*Lx)  ## generation length
L = 10 ## number of diallelic loci
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

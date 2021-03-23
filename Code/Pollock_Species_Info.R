#!/usr/bin/env Rscript
## Dani & Robin
## Population Dynamics
## Genetic and Demographic Independence
#:::::::::::::::::::::::::::::::::::::::::::::::::
# Simulations
#:::::::::::::::::::::::::::::::::::::::::::::::::
# Required Libraries & Functions

#:::::::::::::::::::::::::::::::::::::::::::::::::
## Define Species Parameters
#:::::::::::::::::::::::::::::::::::::::::::::::::
# Pollock Info
# Walleye pollock, Theragra chalcogramma

alpha = 4  # age at maturity, assumed fixed :: 1st mature: 3-4 yr (11), 50%: 4.9 yr/44cm (19)  Appendix Table 1A.2a: Ecological information by life history stage for GOA pollock STOCK ASSESSMENT 2019
MaxAge = 15  ## oldest age for reproduction, after which all die ## Pollock can live longer
####### For initilaistion
# Make survival variable (in function) - SurvivalVariable
# Set up these to initilize pops
Survival = c(1-0.48, 1-0.48, 1-0.48, 1-0.37, rep(1-0.3,c(MaxAge-alpha)))  ## constant survival after age 5 (50% maturity) -> these values are from Table 1.23.
Lx = cumprod(Survival)/Survival  ## cumulative survival
# Number of recruits is variable (in function), but set here to initialize
N1 = 50000
Nx = round(N1*Lx)  ## stable age composition
NTot = sum(Nx)  #sum(Nx), total population of 1+ individuals
AdultN = sum(Nx[alpha:MaxAge])
Fa = c(rep(0, alpha), rep(300, c(MaxAge - alpha))) # Fecundity Oviparous, high fecundity (385- 662⋅103 ) eggs (11) , 1.1-7.2 °C at depth(11)
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
#
bxprime = Bx *adjust
# check
Eggs2 = bxprime*Nx
Teggs2 = sum(Eggs2)  ## eggs required to produce desired number of recruits (N1)
Teggs2*M_a  ## should equal N1

## adjust
MAlpha = M_a # Gompertz alpha
BHBeta = M_a # Becomes 1-B

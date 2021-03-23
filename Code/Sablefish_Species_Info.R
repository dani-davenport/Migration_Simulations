#!/usr/bin/env Rscript
## Dani & Robin
## Population Dynamics
## Genetic and Demographic Indepedance 
#:::::::::::::::::::::::::::::::::::::::::::::::::
# Simulations 
#::::::::::::::::::::::::::::::::::::::::::::::::: 
# Required Libraries & Functions

#:::::::::::::::::::::::::::::::::::::::::::::::::
# # SableFish Info
alpha =   6 # between 5 & 7 for Sablefish age at maturity, assumed fixed
MaxAge =   100 ## oldest age for reproduction, after which all die for sable fish oldest fish caught was 102
####### For initilaistion 
# Make survival variable (in function) - SurvivalVariable
# Set up these to initilize pops
Survival = c(rep(c(1-0.07), 100))  ## differs for male and female in Sablefish constant survival but stable in past assessments, pg45 Stock Assessment, also see (Schirripa, 2007).
Lx = cumprod(Survival)/Survival  ## cumulative survival
# Number of recruits is variable (in function), but set here to initialize 
N1 = 5000 
Nx = round(N1*Lx)  ## stable age composition
NTot = sum(Nx)  #sum(Nx), total population of 1+ individuals
AdultN = sum(Nx[alpha:MaxAge])
GrpAge = 40 # age at which lx drops below 5% of its value at age at first reproduction.
Fa = sapply(list(c(rep(0, alpha), seq(alpha+1, GrpAge), rep(GrpAge, MaxAge-GrpAge))), function(L) 0.000003315*L^3.27264) # weight at age, Table 11, i used female weight exponent a = 3.27264, weight-length coeffiencet = 0.000003315. Figure 51 pg. 162 of the stockassessment shows that Sablfish follow VBGF curve where they dont grow much bigger after the age of 10 years
##  in Sablefish Fecundity eggs/kilogram intercept = 1, slope = 0, that is it appears to be  linearly proportional to weight
Bx = Fa*Nx  ## relative number of offspring produced by each age class
FF = Bx/(sum(Bx))   ## F rescaled to sum to 1 = probability an offspring has parent of age x
units = seq(1:MaxAge)
G = sum(units*Fa*Lx)/sum(Fa*Lx)  ## generation length
L = 10 ## number of diallelic loci

# Popultion Dynamics
# = 0.1 # Gompertz alpha
#BHBeta = 0.1 # Gompertz beta, this is a = 1 - B
##
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

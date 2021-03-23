#!/usr/bin/env Rscript

Reproduce <- function(ParGenos, n_recruits){ ## This is monoecious with mating by age-specific fecundity
  newborns = matrix(0, n_recruits, L)
  # construct matrices of the parental genotypes for each offspring
  # don't want the first col which are ages of parents, don't want last column which is pop origin marker
  newborns =  sapply(1:nrow(newborns), function(x){
    rbinom(n = L, size = 1, ParGenos[[1]][x,2:c(L+1)] / 2) +
      rbinom(n = L, size = 1, ParGenos[[2]][x,2:c(L+1)] / 2) # offspring are a combination of the two parents
  })  #end x
  return(t(newborns))
} #end function


################################
## 1. Intialise popultions
# np = number of popultions  to initialise
# outputs a list length(np), each element conatins a dataframe nrow= NTot ncol = L+3
InitializePops<- function(np){ # function start
  PopsList = list()
  for (j in 1:np){
    ## Initialize population as 100 heterozygotes at L loci
    Genos = matrix(1,NTot,L)
    NSum = cumsum(Nx)  ##running sum of Nx
    Ages = rep(1,NTot)
    for (i in 2:MaxAge) {
      Ages[(NSum[i-1]+1):NSum[i]] = i
    }   ## end for i
    Genos = cbind(Ages,Genos, c(rep(j, nrow(Genos))))  ## add ages for each individual
    PopsList[[j]]<- Genos # Age + Genotypes + Population of Origin
  } ## end for j
  return(PopsList)
}# end InitializePops


# 2. Run Senario
RunSenario_Burnin<- function(PopsList, NYears_Burnin){
  # Burn-in
  BigList = list()
  NYears = NYears_Burnin
  for (j in 1:length(PopsList)){
    BigP = rep(1,L)  ## dummy vector of allele frequencies, year/cycle tracking
    LittleP = rep(1,L) ## dummy vector of allele frequencies, generational tracking
    TotalPopSize = c()
    Genos = PopsList[[j]]

    for (k in 1:NYears) {   ## sampling period
      # Step 1. Total # of Adults
      # realized (Nx) the realized number of individuals alive at time point k
      # this calculates the total number fish in the population/age composition of the population (Nx)
      Nx = seq(1:MaxAge)
      for (i in 2:MaxAge) {
        Nx[i] = sum(Genos[,1]==i)
      } ### end for i
      
      NEggs1 = sum(Bx * Nx) ###  total egg production at stable N
      # mortality alpha
      ## mortality rate required to produce N1 recruits, given nominal bx
      M_a= (N1/NEggs1)
      adjust = (N1/NEggs1)/M_a  ## adjust > 1 indicates bx is too low, given M
      bxprime = Bx *adjust
      # check
      Eggs2 = bxprime*Nx
      Teggs2 = sum(Eggs2)  ## eggs required to produce desired number of recruits (N1)
      
      # final Eggs
      
      NEggs2 = sum(bxprime * Nx)

      # Step 3. Determine N1 using Gompertz (Okamoto Paper, Eqn 3) with environment effect on recruitment
      n_recruits = round(MAlpha*NEggs2^(1-BHBeta)*exp(rnorm(1, 0, 0.25))) # natural log, environment effect of recruitment
      print(n_recruits)
      
      ## create newborn offspring
      ## 1. Get Parents, a vector of indicies
      M = cbind(c(1:length(Genos[,1])),sapply(Genos[,1], function(x) bxprime[x])) # weighted by age specific fecundity
      # Rescale weights
      M[,2] <- M[,2]/sum(M[,2])
      # Parental Genotypes
      ParGenos = list(Genos[sample(M[,1], size = n_recruits, replace = TRUE, prob = M[,2]),],
                      Genos[sample(M[,1], size = n_recruits, replace = TRUE, prob = M[,2]),]) # list of two lots of parents
      # Make recruits
      recruits = Reproduce(ParGenos, n_recruits)

      HZ = colMeans(recruits[,1:L])/2

      if (k==1) {
        BigP = HZ } else {
          BigP = rbind(BigP,HZ)  }  ## end if

      newAge = rep(1,nrow(recruits))
      recruits = cbind(newAge,recruits, c(rep(j, nrow(recruits)))) # age, genotype, population of origin


      ## Delete individuals in oldest age class
      Genos = subset(Genos, Genos[,1]<MaxAge)

      ## Everybody ages 1 year
      Genos[,1] = Genos[,1]+1

      ## Some die
      #NSurvive = NTot - N1
      NSurvive = which(rbinom(nrow(Genos), 1, Survival[1]) == 0) # does not account for age
      Genos = Genos[-NSurvive,]

      ## add new recruits
      Genos <- rbind(recruits,Genos)

      LittleP = rbind(LittleP, colMeans(Genos[,2:c(L+1)])/2)

      PopsList[[j]] <- Genos # the number of recruits varies, all "recruits" or "newborns" survive
      # Track popultion size
      TotalPopSize = c(TotalPopSize, nrow(Genos)) # includes new recruits
    }  ## end for k
    BigList[[j]]<- list(Genotypes=PopsList[[j]], AlleleFreqs=BigP, TotalPopSize = TotalPopSize, AlleleFreqsPop = LittleP)
  } ## end for j
  return(BigList)
} # end function

RunSenario_Migrate<- function(PopsList, lambda = lambda, NYears_Migrate) {
  #### MIGRATE
  print("start migration")
  #PopsList = PopsList
  NYears = NYears_Migrate

  for (t in 1:NYears){
    Recruits_List = list()
    for (j in 1:length(PopsList)){
      Genos = PopsList[[j]]$Genotypes
      # Step 1. Total # of Adults
      Nx = seq(1:MaxAge)
      for (i in 2:MaxAge) {
        Nx[i] = sum(Genos[,1]==i)
      } ### end for i
      # Step 1. Fix eggs
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
      #Step 2. 
      NEggs2 = sum(bxprime * Nx)  ## eggs required to produce desired number of recruits (N1)

      # Step 3. Determine N1 using Gompertz (Okamoto Paper, Eqn 3) with environment effect on recruitment
      n_recruits = round(MAlpha * NEggs2^(1-BHBeta) * exp(rnorm(1, 0, 0.25))) # natural log, environment effect of recruitment

      
      ## create newborn offspring
      ## 1. Get Parents, a vector of indicies
      ## 1. Get Parents, a vector of indicies
      M = cbind(c(1:length(Genos[,1])),sapply(Genos[,1], function(x) bxprime[x])) # weighted by age specific fecundity
      # Rescale weights
      M[,2] <- M[,2]/sum(M[,2])
      # Parental Genotypes
      ParGenos = list(Genos[sample(M[,1], size = n_recruits, replace = TRUE, prob = M[,2]),],
                      Genos[sample(M[,1], size = n_recruits, replace = TRUE, prob = M[,2]),]) # list of two lots of parents
      recruits = Reproduce(ParGenos, n_recruits)
      Recruits_List[[j]] <-  cbind(recruits, rep(j, nrow(recruits))) # track population of origin
    } # end j

    # Initiate Migration
    new_recruits = Migrate(Recruits_List, lambda) # returns a list with genotype and chromosome information

    for (i in 1:length(PopsList)){
      Genos <- PopsList[[i]]$Genotypes

      LittleP <- PopsList[[i]]$AlleleFreqsPop
      BigP <- PopsList[[i]]$AlleleFreqs # this will join it to existing AF
      livefish <- new_recruits[[i]][,1:L]
      BigP <- rbind(BigP, colMeans(livefish)/2)
      # New recuits age 1
      newAge = rep(1, nrow(new_recruits[[i]]))
      recruits = cbind(newAge,new_recruits[[i]])  # ages, genotypes and chromosome information joined

      ## Delete individuals in oldest age class
      Genos = subset(Genos, Genos[,1]<MaxAge)
      ## Everybody ages 1 year
      Genos[,1] = Genos[,1]+1
      ## Some die

      NSurvive = which(rbinom(nrow(Genos), 1, Survival[1]) == 0) # does not account for age
      Genos = Genos[-NSurvive,]

      ## set colnames or datasets wont join
      #recruits = setNames(data.frame(recruits), c(colnames(Genos)))
      ## Add new recruits
      Genos = rbind(recruits,Genos)

      # track proportion of alleles in each pop that are from migrant populations
      #migrantAlleles = getPropMigrantAllele(Genos, i)
      # track population of origin
      migrantPop = getPropMigrants(Genos, i)

      PopsList[[i]]$Genotypes <-  Genos


      LittleP = rbind(LittleP, colMeans(Genos[,2:c(L+1)])/2) # population Allele Frequency

      PopsList[[i]]$AlleleFreqs <- BigP
      PopsList[[i]]$AlleleFreqsPop <- LittleP
      PopsList[[i]]$TotalPopSize <- c(PopsList[[i]]$TotalPopSize, nrow(Genos)) # includes new recruits
      #PopsList[[i]]$AlleleExchange <- c(PopsList[[i]]$AlleleExchange, migrantAlleles)
      PopsList[[i]]$PopOrigin <- c(PopsList[[i]]$PopOrigin, migrantPop) # last column
    } # end i
  } # end t
  return(PopsList)
}


################################
# Migration between Pops
# newborns will switch between popultions after they are created
Migrate = function(recruits, lambda){
  # The function takes (1) recruits, a list of lists with [[]] a df of genotypes of recruits
  #                    (2) r, a migration rate (both pops with the same migration rate)
  mrecruits_pop1 = which(rbinom(c(1:nrow(recruits[[1]])), 1, lambda) == 1) #  migration is a bionomial draw
  mrecruits_pop2 = which(rbinom(c(1:nrow(recruits[[2]])), 1, lambda) == 1) #  migration is a bionomial draw
  # forcing a % of migration, per generation
  #mrecruits_pop1 = sample(recruits[[1]], nrow(recruits[[1]])*lambda, replace = FALSE) #prob =
  #mrecruits_pop2 = sample(recruits[[2]], nrow(recruits[[2]])*lambda, replace = FALSE)
  new_recruits_pop1 = rbind(recruits[[2]][mrecruits_pop2,], recruits[[1]][-mrecruits_pop1,])
  new_recruits_pop2 = rbind(recruits[[2]][-mrecruits_pop2,], recruits[[1]][mrecruits_pop1,])
  #new_recruits_pop1 = merge(recruits[[2]][mrecruits_pop2,],recruits[[1]][-mrecruits_pop1,])
  #new_recruits_pop2 = merge(recruits[[2]][-mrecruits_pop2,], recruits[[1]][mrecruits_pop1,])
  return(list(new_recruits_pop1, new_recruits_pop2))

} # end function Migrate

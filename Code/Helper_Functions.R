#!/usr/bin/env Rscript
#:::::::::::::::::::::::::::::::::::::::::::::::::
# Metrics
#:::::::::::::::::::::::::::::::::::::::::::::::::



normalize <- function(x){
  return((x-min(x)) / (max(x)-min(x)))}

# get Fst function (per loci, over time)
GetFst<- function(v, N){ #function(v, N) v is a matrix of af at one time point, and the comparisions to make, N is a vector of population sizes
  FstVec= c()
  # Ny and Nx, size of the popultions at time(p) ## need to track pop size...
  for (l in 1:ncol(v)){ # # of comparisions
    x = as.numeric(as.character(v[1,l]))
    y = as.numeric(as.character(v[2,l]))
    #Nx = N[1,] # pop size 1 at time t
    #Ny = N[2,] # pop size 2 at time t
    ### Wrights Varience Fst 
    meanAF = (x+y)/2 
    varp = (((x-meanAF)^2) + ((y-meanAF)^2)) # population variance Var[p]
    phat = (x+y)/2 
    Fst = varp/phat*(1-phat)
    FstVec[l] = Fst
  } # end l
  return(FstVec)
  #return(Fst)
} # end function 


# Get the fraction of migratnts in the recipent popultion 
## migration goes both ways 
## get the proportion of # 2 in pop 1
## get the proportion of # 1 in pop 2 
## get for the whole pop (recruits and adults)

getPropMigrants<- function(Genos, pop){ 
  # function returns the proportion of migrants as a function of the total popultion size at time 
  totalpopsize = nrow(Genos)
  if (pop ==1) {
    prop = length(which(Genos[,L+2]==2))/totalpopsize # Genos has ages in first row + L loci
  }else{
    prop = length(which(Genos[,L+2]==1))/totalpopsize  
  } # end if else
  return(prop)
} # end getPropMigrant


slideFunct_Mean <- function(vec, window, step){
  varNSlide = c()
  indx = as.numeric(seq(from=1, to=length(vec), by=step))
  for (i in as.numeric(indx)){
    #print(i)
    #print(i+window-1)
    #meanNSlide <- c(meanNSlide, mean(vec[i:c(i + window - 1)]))
    varNSlide <- c(varNSlide, var(vec[i:c(i + window-1)]))
    
  } # end i 
  #return(meanNSlide)
  return(varNSlide)
} # end function 
## Min MaxN in sliding window 
minMaxN = c()
slideFunct_MinMax<- function(data, window, step){
  indx = as.numeric(seq(from=1, to=(length(data)), by=step))
  print(indx)
  for (i in as.numeric(indx)){
    minMaxN <- rbind(minMaxN, c(min(data[i:c(i + window - 1)]), max(data[i:c(i + window - 1)])))
  } # end i 
  return(minMaxN)
} # end function 


getDiff<- function(x) abs(diff(x)) #/x[-length(x)] # differences in N overtime, AF overtime 
# parametric F, change in AF between each time point, where P is the frequeny of the allele in generation t
getDiffAF<- function(x) (diff(x)^2)/(x[-length(x)]*(1-x[-length(x)]))

GetFst<- function(v, N){ #function(v, N) v is a matrix of AFS at one time point, and the comparisions to make, N is a vector of population sizes
  FstVec= c()
  # Ny and Nx, size of the popultions at time(p) ## need to track pop size...
  for (l in 1:ncol(v)){ # # of comparisions
    x = as.numeric(as.character(v[1,l]))
    y = as.numeric(as.character(v[2,l]))
    #Nx = N[1,] # pop size 1 at time t
    #Ny = N[2,] # pop size 2 at time t
    ### Wrights Varience Fst 
    meanAF = (x+y)/2 
    varp = (((x-meanAF)^2) + ((y-meanAF)^2)) # population variance Var[p]
    phat = (x+y)/2 
    Fst = varp/phat*(1-phat)
    FstVec[l] = Fst
  } # end l
  return(FstVec)
  #return(Fst)
} # end function 

# gets the data for Fst to plot
SamplePopsFS<- function(data_in, numReps, j){
  Fst_Time = c()
  r = numReps
  #print(subSample)
  #for (i in subSample){
    #print(i)
    data = data_in[[j]]
  # get each rep

    Ntime = nrow(data[[1]]$AlleleFreqsPop)
  # tracking Fst across time
  # two pops only 
  for (t in 1:Ntime)  {
    pairwiseAF= c()
    Nsize = c() # list of pop sizes at time t
    for (p in 1:length(data)){
      pairwiseAF = rbind(pairwiseAF, c(data[[p]]$AlleleFreqsPop[t,])) # loci cols, pops in rows
      Nsize = c(Nsize, data[[p]]$TotalPopSize[t])
    } # end p 
    N = combn(Nsize, 2)
    Fst = GetFst(pairwiseAF, N)
    # add a column with "rep" number for plotting 
    Fst = c(t,Fst, c(as.character(r)))
    ##} # end l 
    Fst_Time = rbind(Fst_Time, Fst) # fst at time t, pop comparisons in rows 
  } # end t
  #} # end i 
  colnames(Fst_Time) = c("Time", paste0("Loci", 1:L), "Rep")
  return(Fst_Time)
} # end function


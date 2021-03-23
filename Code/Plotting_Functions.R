## Code to make Plots 
pkgs_plots = c("dplyr", "ggplot2", "tidyr", "ggpubr","purrr", "boot")
lapply(pkgs_plots, library, character.only = TRUE)  
source("Plotting_Functions")
source("Code/Helper_Functions.R")
# plot colors
under.the.sea.palette <- c("#1C2344", "#163D7D", "#118386", "#89D6D6", "#6E783D", "#FC564F", "#901546", "#FFD521", "#010101")

# you can subset the number of reps to test the effects on the results
numRep_List<- c(800) #, 400, 200, 100, 50, 25)

# read in the R data 
# make empty list to take data 
input <- list()
# get list of files containing data
list.files(pattern = "\\.dbf$")
for(i in 1:length(list.files)){
  load(list.files[[i]])
  input[[i]] <- Herring
  }



# code to take a random subset of sims
# output is a list of length 3, where each element represents reps from either 0.05, 0.1, 0.2 migration rate
# input_sample<- lapply(Hake, function(d) lapply(numRep_List, function(x) sample(d, x)))
# input =  input_sample# output from sims
  
# will work on any length of numRep_List
# Make Plots 
lapply(1:length(lambdaList), function(x) 
  makeSimPlots(x, data = input[[x]], numRepList= numRep_List))

# input = input_sample

# make correlation N plots (currently only works for a single rep list e.g. numRep_List = c(1000))
make_corrN_before(input, lambdaList) #logN
make_corrN_after(input, lambdaList)  #logN
# scatter
make_corrN_before_after_scatter(input, lambdaList, numRepList= numRep_List)
# 
# ??# still to do?
# # make correlation AF plots (currently only works for a single rep list e.g. numRep_List = c(1000))
# make_corrAF_before(input, lambdaList) ##
# make_corrAF_after(input, lambdaList) ##
# # scatter
# make_corrAF_before_after_scatter() ##
N = 10
subSample = sample(1:numRep_List[[1]], N, replace = F) # Need a representative Sample, this makes an index which will be used to make plots of a subsample of populaiton replictaes for Fst and P1 v P2 N plots

getN_Normalised_P1vP2(input, N = 10, numRepList = numRep_List, lambdaList, subSample = subSample) # only uses first item in the list as input
# make Fst plots, uses subsample so the plots are from the same "replicates" as used in the getN_Normalised plot                                  
makeFstPlot(input, N = 10, numRepList = numRep_List, lambdaList, subSample = subSample)                                 


#lapply(1:length(lambdaList), function(x) 
#  makeSimPlotsMigrantAlleles(x, data = Hake_Test[[x]], numRepList= numRep_List))


#:::::::::::::::::::::::::::::::::::::::::::::::::
################### PLOT ########################

total_cycles = BurnIn_Reps + Collect_Reps + Migrate_Reps

# Plotting the results  between replicate runs of the above senario

changeN_List<- list()
coorNBefore<- list()
coorNAfter<- list() # These will only work for one numRepList, not multiple different numbers of replicates 
coorAFBefore<- list()
coorAFAfter<- list()

makeSimPlots<- function(migRate, dataList, numRepList){
  migrationStart = Collect_Reps + BurnIn_Reps + 100 # some period of time after migration to start collecting data
  migrationEnd = total_cycles
  
  ## numReps is a list 
  ## data is a list 
  ## migRate is an index of migration rate, used for labelling and indexing
  migRate = lambdaList[migRate]
  for (i in 1:length(numRepList)) {

    
    data = dataList[[i]]
    numRep = numRepList[i]
    # 1. Get changes of N (first differences) over time, across replicates 
    changeN=  data.frame()
    # where change in N is averaged across replicates using getDiff
    changeN_Mig= data.frame()
    for (r in 1:numRep){
      diffN= sapply(list(data[[r]][[1]]$TotalPopSize, data[[r]][[2]]$TotalPopSize), function(x) getDiff(x))
      diffN= as.data.frame(cbind(rep(r, nrow(diffN)), seq(1:nrow(diffN)), diffN))
      colnames(diffN)<- c("Rep", "Time", "Pop1", "Pop2")
      changeN_Mig = rbind(changeN_Mig, diffN)
    } # end r 
    ## plot 
    changeN_Mig %>%
      select(Rep, Time, Pop1, Pop2) %>%
      gather(key = "Pop", value = "Change", c(-Time, -Rep)) %>%
      group_by(Time, Pop) %>% # all loci avergaged
      summarize(mean_change = mean(Change, na.rm = TRUE)) %>%  # the average change in N 
      ggplot(., aes(x = Time, y = mean_change)) +  #, linetype = factor(Pop)
      geom_line(aes(color = factor(Pop))) +
      #aes(group = factor(Pop)),
      #          size=0.8) + 
      theme_bw(base_size = 18) + 
      labs(y = "Mean Absolute Change N") +
      scale_color_manual(values = c("grey3","grey67")) + 
      guides(color = "none")
    #scale_x_continuous() 
    #facet_wrap(~Rep, ncol =1) + guides(fill = FALSE, color = FALSE) 
    ggsave(paste0(paste(numRep, migRate, sep = "_"), "_Avg-Change-N.png"),width = 30, height = 20, units = "cm")
    
    #2.    Mean, Varience in N, sliding window 
    # Get mean in window size, set to generation length 
    window = 10
    slide = 10
    varNslideDF = data.frame()
    for (r in 1:numRep){
      out = sapply(list(data[[r]][[1]]$TotalPopSize, data[[r]][[2]]$TotalPopSize), function(x) slideFunct_Mean(x, window, slide)) 
      out = as.data.frame(cbind(rep(r, nrow(out)), seq(1:nrow(out)), out))
      colnames(out)<- c("Rep", "Time", "Pop1", "Pop2")
      varNslideDF  = rbind(varNslideDF , out)
    } # end r 
    # Plot the varience in a sliding window, avg across replicates 
    
    varNslideDF %>%
      select(Rep, Time, Pop1, Pop2) %>%
      gather(key = "Pop", value = "Var_N_Slide", c(-Time, -Rep)) %>%
      group_by(Time, Pop) %>% 
      summarize(mean_change_var = mean(Var_N_Slide, na.rm = TRUE)) %>%  
      ggplot(., aes(x = Time, y = mean_change_var, linetype = factor(Pop))) + 
      geom_line(aes(group = factor(Pop)),
                size=0.8) + #aes(color = Pop)) + 
      theme_bw(base_size = 18) + 
      labs(y = "Mean Variance N") + 
      scale_color_manual(values = c("grey3","grey67")) +
      guides(color = "none")
    #scale_x_continuous(limits=c(0, 75), breaks = c(1,25,75)) 
    #ggsci::scale_color_jco() #+ 
    #scale_color_manual(values = c("darkred", "steelblue")) + 
    #facet_wrap(~Rep, ncol =1) + guides(fill = FALSE, color = FALSE) 
    ggsave(paste0(paste(numRep, migRate, sep = "_"),"_Avg-Var-N-SlidingWindow-10.png"),width = 30, height = 20, units = "cm")
    # # get min max N 
    # minMaxNslideDF = data.frame()
    # window = 10
    # slide = 10
    # for (r in 1:numRep){
    #   out = sapply(list(data[[r]][[1]]$TotalPopSize, data[[r]][[2]]$TotalPopSize), function(x) slideFunct_MinMax(x, window, step)) 
    #   print(head(out))
    #   colnames(out)<- c("Rep", "Time", "Pop1-Min", "Pop1-Max", "Pop1-Min", "Pop1-Max")
    #   minMaxNslideDF  = rbind(varNslideDF , out)
    # } # end r 
    
    #:::::::::::::::::::::::::::::::::::::::::::::::::
    ## Correlations
    #:::::::::::::::::::::::::::::::::::::::::::::::::
    # 3A. The correlation between popultation size N between Pop1 and Pop2, the log10 of N, then calculate the correlation 
    # coorsN<- c()
    # ### before migration 
    # for (r in 1:numRep){
    #   xN= sapply(list(data[[r]][[1]]$TotalPopSize[BurnIn_Reps:Collect_Reps+BurnIn_Reps], data[[r]][[2]]$TotalPopSize[BurnIn_Reps:Collect_Reps+BurnIn_Reps]), function(x) x)
    #   xN= as.data.frame(cbind(rep(r, nrow(xN)), seq(1:nrow(xN)), xN))
    #   colnames(xN)<- c("Rep", "Time", "Pop1", "Pop2")
    #   xCor = log10(xN$Pop1)
    #   yCor = log10(xN$Pop2)
    #   coorsN <- c(coorsN,cor(xCor,yCor)) 
    # } # end r 
    #  
    # # Plot Coors 
    # ggplot2::qplot(coorsN, geom="histogram", binwidth = 0.05) +  theme_bw() + 
    #   theme(axis.text=element_text(size=20)) + 
    #   labs(x = "Pearson r", fill = "lambda") + 
    #   scale_fill_manual(values = under.the.sea.palette) + 
    #   scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
    #   scale_y_continuous(expand = c(0, 0)) 
    # 
    # ggsave(paste0(paste(numRep, migRate, sep = "_"),"_Correlation-log10N-Before-Migration.png"), width=30, height=20, units="cm")
    # 
    #3B. Corelation in N after migration 
    # coorsN<- c()
    # 
    # for (r in 1:numRep){
    #   xN= sapply(list(data[[r]][[1]]$TotalPopSize[migrationStart:migrationEnd], data[[r]][[2]]$TotalPopSize[migrationStart:migrationEnd]), function(x) x)
    #   xN= as.data.frame(cbind(rep(r, nrow(xN)), seq(1:nrow(xN)), xN))
    #   colnames(xN)<- c("Rep", "Time", "Pop1", "Pop2")
    #   x = xN$Pop1
    #   y = xN$Pop2
    #   coorsN <- c(coorsN,cor(x,y)) 
    # } # end r 
    # coorNAfter[[i]]<- coorsN 
    # ggplot2::qplot(coorsN, geom="histogram", binwidth = 0.05) +   theme_bw() + 
    #   theme(axis.text=element_text(size=20)) + 
    #   labs(x = "Pearson r", fill = "lambda") + 
    #   scale_fill_manual(values = under.the.sea.palette) + 
    #   scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
    #   scale_y_continuous(expand = c(0, 0)) 
    # 
    # ggsave(paste0(paste(numRep, migRate, sep = "_"),"_Correlation-N-After-Migration.png"), width=30, height=20, units="cm")
    # 
    # #:::::::::::::::::::::::::::::::::::::::::::::::::
    # 4. Allele Frequencies 
    # Get changes of AF (parameteric F) over time, across replicates 
    # All replicates,  this gets the average change in AF across replicates at each time point. 
    changeAF <- data.frame()
    for (r in 1:numRep){
      for ( l in 1:L){ # for each loci 
        diffAF= sapply(list(unname(unlist(data[[r]][[1]]$AlleleFreqsPop[, l])), unname(unlist(data[[r]][[2]]$AlleleFreqsPop[, l]))), function(x) getDiffAF(x))
        diffAF= as.data.frame(cbind(as.numeric(rep(l, nrow(diffAF))), rep(r, nrow(diffAF)),seq(1:nrow(diffAF)), diffAF))
        colnames(diffAF)<- c("Loci", "Rep", "Time", "Pop1", "Pop2")
        changeAF <- rbind(changeAF, diffAF)
        #changeAF[r,l]<- c(cor(xDif , yDif))
      } # end l
    } # end r
    
    changeAF %>%
      select(Loci, Rep, Time, Pop1, Pop2) %>%
      #filter(Rep == 1) %>%
      gather(key = "Pop", value = "Change", c(-Time, -Rep, -Loci)) %>%
      mutate(Loci = factor(Loci, levels=c(seq(1:L))))%>%
      group_by(Time, Pop) %>% # all loci avergaged
      summarize(mean_change = mean(Change, na.rm = TRUE)) %>%  # the average change in AF 
      ggplot(., aes(x = Time, y = mean_change, linetype = factor(Pop))) + 
      geom_line(aes(group = factor(Pop)),
                size=0.8) + 
      theme_bw(base_size = 18) + 
      labs(x="Time", y = "F") +
      scale_color_manual(values = c("grey3","grey67")) + 
      guides(color = "none")
    #scale_x_continuous(breaks = c(0,50,250,750), limits=c(0, 750)) 
    ggsave(paste0(paste(numRep, migRate, sep = "_"),"_Avg-Change-AF-Population.png"),width = 30, height = 20, units = "cm")
    
    #5A. Get correlation in AF, before migration
    coorsAF<- data.frame()
    for (r in 1:numRep){
      for (l in 1:L){ # for each loci 
        xAF= sapply(list(unname(unlist(data[[r]][[1]]$AlleleFreqsPop[BurnIn_Reps:Collect_Reps+BurnIn_Reps, l])), unname(unlist(data[[r]][[2]]$AlleleFreqsPop[BurnIn_Reps:Collect_Reps+BurnIn_Reps, l]))), function(x) x)
        xAF= as.data.frame(cbind(rep(l, nrow(xAF)), rep(r, nrow(xAF)), seq(1:nrow(xAF)), xAF))
        colnames(xAF)<- c("Loci", "Rep", "Time", "Pop1", "Pop2")
        x = xAF$Pop1
        y = xAF$Pop2
        coorsAF[r,l]<- cor(x,y) } # end l
    } # end r 
    coorAFBefore[[i]]<- coorsAF
    # Plot Coors 
    # Robin says do correlation bewteen each locus in each replicate 
    # take the average correlation across loci at each time point 
    #coorsAFAvg = rowMeans(coorsAF)
    coorsAF_loci = unlist(coorsAF)
    qplot(coorsAF_loci, geom="histogram", binwidth = 0.2) +  theme_bw() + 
      theme(axis.text=element_text(size=20)) + 
      labs(x = "Pearson r", fill = "lambda") + 
      scale_fill_manual(values = under.the.sea.palette) + 
      scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0)) 
    ggsave(paste0(paste(numRep, migRate, sep = "_"),"_Correlation-AF-Loci-Before-Migration.png"), width=30, height=20, units="cm")
    
    #5B. Correlation in AF over time, after migration 
    coorsAF<- data.frame()
    for (r in 1:numRep){
      for (l in 1:L){ # for each loci 
        xAF= sapply(list(unname(unlist(data[[r]][[1]]$AlleleFreqsPop[migrationStart:migrationEnd, l])), unname(unlist(data[[r]][[2]]$AlleleFreqsPop[migrationStart:migrationEnd, l]))), function(x) x)
        xAF= as.data.frame(cbind(rep(l, nrow(xAF)), rep(r, nrow(xAF)), seq(1:nrow(xAF)), xAF))
        colnames(xAF)<- c("Loci", "Rep", "Time", "Pop1", "Pop2")
        x = xAF$Pop1
        y = xAF$Pop2
        coorsAF[r,l]<- cor(x,y) } # end l
    } # end r 
    coorAFAfter[[i]]<- coorsAF
    # Plot Coors 
    # take the average correlation across loci at each time point 
    #coorsAFAvg = rowMeans(coorsAF)
    coorsAF_loci= unlist(coorsAF)
    qplot(coorsAF_loci, geom="histogram", binwidth = 0.2) +  theme_bw() +  theme(axis.text=element_text(size=20))
    ggsave(paste0(paste(numRep, migRate, sep = "_"),"_Correlation-AF-After_Migration-Loci.png"), width=30, height=20, units="cm")
    
    #:::::::::::::::::::::::::::::::::::::::::::::::::
    # 6. Plot avg. fraction of  migrant genes across reps
    df_AE= data.frame()
    
    for (r in 1:numRep){
      
      AlleleExchangeDF = cbind(Pop_1 = data[[r]][[1]]$PopOrigin, Pop_2 = data[[r]][[2]]$PopOrigin, Time = seq(1:length(data[[r]][[1]]$PopOrigin))) # where time is the number of cycles when migration occoured 
      AElong = as.data.frame(AlleleExchangeDF) %>% tidyr::pivot_longer(cols = starts_with("Pop"), 
                                                                       names_to = "Pop", 
                                                                       names_prefix = "Pop_", 
                                                                       values_to = "AE", 
                                                                       values_drop_na = F)
      if (r == 1){
        df_AE= AElong}else{
          df_AE = cbind(df_AE, AElong$AE)
        }# end if
    } # end r
    
    names(df_AE) <- c("Time", "Pop", paste0("Rep_", seq(1:numRep)))
    df_AE %>%
      #select(-Time, -Pop) %>% 
      mutate(MEAN = rowMeans(.[,-(1:2)])) %>%
      # change time for x-axis
      ggplot(data = ., aes(x=Time, y= MEAN, group = Pop, color=Pop)) +
      geom_line() + 
      #geom_point(size = 0.5) + 
      #scale_color_brewer(palette="Paired") + 
      labs(x="Time", y = "Propotion Migrant per Population") +
      #scale_x_discrete(limits=c(0,500,550)) +
      theme_bw(base_size = 18) + 
      scale_x_continuous(breaks = c(seq(0, Migrate_Reps, by = 100)), 
                         labels = c("Start", seq(100,Migrate_Reps, by = 100))) + 
      scale_color_manual(values = c("grey3","grey67")) + 
      ylim(0,0.50) +
      guides(fill = FALSE, color = FALSE)
    
    ggsave(paste0(paste(numRep, migRate, sep = "_"), "_Avg_Proportion_of_Migrant_Alleles-After-Migration.png"), width=30, height=20, units="cm")
    
  } # end i
  return(NULL)
} # end function


make_corrN_before_after_scatter = function(input, lambdaList, numRepList= numRep_List){
  numRep = numRep_List[[1]]
  coorN_before = lapply(1:length(lambdaList), function(x) 
    coorN1(x, data = input[[x]], numRepList= numRep_List))
  coorN_after = lapply(1:length(lambdaList), function(x) 
    coorN2(x, data = input[[x]], numRepList= numRep_List))
  # after
  y = as.data.frame(cbind(1:numRep, do.call(cbind, lapply(coorN_after, function(x) do.call(cbind, x)))))
  colnames(y) <- c("NumRep", lambdaList)
  y = tidyr::pivot_longer(data = y, cols = -NumRep, names_to = "MigRate", values_to = "Corr")
  # before
  x = as.data.frame(cbind(1:numRep, do.call(cbind, lapply(coorN_before, function(x) do.call(cbind, x)))))
  colnames(x) <- c("NumRep", lambdaList)
  x = tidyr::pivot_longer(data = x, cols = -NumRep, names_to = "MigRate", values_to = "Corr")
  out = lapply(list(x,y), function(x) x %>%  
                 group_by(MigRate) %>%
                 nest() %>% 
                 mutate(boot_res = map(data,
                                       ~ boot(data = .$Corr,
                                              statistic = function(x, i) median(x[i]),
                                              R = 1000)),
                        boot_res_ci = map(boot_res, boot.ci, type = "perc"),
                        median = map(boot_res_ci, ~ .$t0),
                        lower_ci = map(boot_res_ci, ~ .$percent[[4]]),
                        upper_ci = map(boot_res_ci, ~ .$percent[[5]]), # 95% CIs
                        n =  map(data, nrow)) %>% 
                 select(-data, -boot_res, -boot_res_ci) %>% 
                 unnest(cols = c(n, median, lower_ci, upper_ci)) %>% 
                 ungroup())
  
  out[[1]]<- tibble::add_column(out[[1]], Time = c(rep("Before", length(lambdaList))))
  out[[2]]<- tibble::add_column(out[[2]], Time = c(rep("After", length(lambdaList))))
  
  
  # make plot
  pd = position_dodge(-0.2)
  p =  out %>%
    purrr::reduce(.,  dplyr::bind_rows) %>%
    ggplot(.,aes(x = MigRate, y= median, color = Time)) + 
    geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=.1, position=pd) +
    geom_point(position=pd) + 
    labs(x = "lambda", y = "median Pearson r (95%CIs)") + 
    theme_bw(base_size = 20) + 
    theme(legend.title=element_blank()) + 
    scale_color_manual(values = c("#010101","#118386"))
  ggsave(paste0(numRep, "_Correlation-medianlog10N-Scatter95CIs.png"),width = 30, height = 20, units = "cm")
}





make_corrN_after = function(input, lambdaList){
  
  coorN_after = lapply(1:length(lambdaList), function(x) 
    coorN2(x, data = input[[x]], numRepList= numRep_List))
  
  
  x = as.data.frame(cbind(1:numRep, do.call(cbind, lapply(coorN_after, function(x) do.call(cbind, x)))))
  colnames(x) <- c("NumRep", lambdaList)
  x = tidyr::pivot_longer(data = x, cols = -NumRep, names_to = "MigRate", values_to = "Corr")
  ggplot2::ggplot(x, aes(Corr, fill = MigRate)) +
    geom_histogram(binwidth = 0.1,  position="dodge") + #, alpha=0.5 position = "identity"
    theme_bw() + 
    theme(text=element_text(size=20)) +
    labs(x = "Pearson r (logN)", fill = "Migration Rate") + 
    scale_fill_manual(values = under.the.sea.palette) + 
    scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  ggsave(paste0(numRep, "_Correlation-log10N-After-Migration.png"),width = 30, height = 20, units = "cm")
  
} # end function


make_corrN_before = function(input, lambdaList){
  
  coorN_before = lapply(1:length(lambdaList), function(x) 
    coorN1(x, data = input[[x]], numRepList= numRep_List))
  
  x = as.data.frame(cbind(1:numRep, do.call(cbind, lapply(coorN_before, function(x) do.call(cbind, x)))))
  colnames(x) <- c("NumRep", lambdaList)
  x = tidyr::pivot_longer(data = x, cols = -NumRep, names_to = "MigRate", values_to = "Corr")
  ggplot2::ggplot(x, aes(Corr, fill = MigRate)) +
    geom_histogram(binwidth = 0.1,  position="dodge") + #, alpha=0.5 position = "identity"
    theme_bw() + 
    theme(text=element_text(size=20)) +
    labs(x = "Pearson r (logN)", fill = "Migration Rate") + 
    scale_fill_manual(values = under.the.sea.palette) + 
    scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  ggsave(paste0(numRep, "_Correlation-log10N-Before-Migration.png"),width = 30, height = 20, units = "cm")
  
} # end function

######## --------
coorN2<- function(migRate, dataList, numRepList){
  migrationStart = Collect_Reps + BurnIn_Reps + 100 # some period of time after migration to start collecting data
  migrationEnd = total_cycles
  
  ## numReps is a list 
  ## data is a list 
  ## migRate is an index of migration rate, used for labelling and indexing
  migRate = lambdaList[migRate]
  for (i in 1:length(numRepList)) {
    data = dataList[[i]]
    numRep = numRepList[i]
    coorNAfter = list()
    coorsN<- c()
    ### before migration 
    for (r in 1:numRep){
      xN= sapply(list(data[[r]][[1]]$TotalPopSize[migrationStart:migrationEnd], data[[r]][[2]]$TotalPopSize[migrationStart:migrationEnd]), function(x) x)
      xN= as.data.frame(cbind(rep(r, nrow(xN)), seq(1:nrow(xN)), xN))
      colnames(xN)<- c("Rep", "Time", "Pop1", "Pop2")
      xCor = log10(xN$Pop1)
      yCor = log10(xN$Pop2)
      coorsN <- c(coorsN,cor(xCor,yCor)) 
      coorNAfter[[i]] <- coorsN
    } # end r 
  } # end i
  return(coorNAfter)
} # end function


######## --------
coorN1<- function(migRate, dataList, numRepList){
  
  ## numReps is a list 
  ## data is a list 
  ## migRate is an index of migration rate, used for labelling and indexing
  migRate = lambdaList[migRate]
  for (i in 1:length(numRepList)) {
    data = dataList[[i]]
    numRep = numRepList[i]
    coorNBefore = list()
    coorsN<- c()
    ### before migration 
    for (r in 1:numRep){
      xN= sapply(list(data[[r]][[1]]$TotalPopSize[BurnIn_Reps:Collect_Reps+BurnIn_Reps], data[[r]][[2]]$TotalPopSize[BurnIn_Reps:Collect_Reps+BurnIn_Reps]), function(x) x)
      xN= as.data.frame(cbind(rep(r, nrow(xN)), seq(1:nrow(xN)), xN))
      colnames(xN)<- c("Rep", "Time", "Pop1", "Pop2")
      xCor = log10(xN$Pop1)
      yCor = log10(xN$Pop2)
      coorsN <- c(coorsN,cor(xCor,yCor)) 
      coorNBefore[[i]] <- coorsN
    } # end r 
  } # end i
  return(coorNBefore)
} # end function



# get abundance for plotting 
## migRate is an index of migration rate, used for labelling and indexing
getN_Normalised_P1vP2 = function(input, N, numRepList, lambdaList, subSample){ # N is the number (index) of replicates you'd like to plot
  # if (time == "BEFORE"){}
  START = BurnIn_Reps + 1
  END = total_cycles
  numRep = subSample
  
  pltN = lapply(1:length(lambdaList), function(x) {
    plotNV(data_in = input[[x]], numReps = numRep) # i normalise the data
  })
  # make the plots
  names(pltN)<- lambdaList
  # collapse the list 
  plist = lapply(1:length(lambdaList), function(i){ 
    x = pltN[[i]]
    name = names(pltN[i])
    x %>%
    tibble::add_column(., MigRate = rep(name, nrow(x))) %>%
    select(-LogN) %>%
    mutate(., Pops =  factor(Pops)) %>%
    mutate(., Rep =  factor(Rep)) %>%
    spread(., key = Pops, value = NormLogN)
  })
    
   plistdf = as.data.frame(do.call(rbind, plist))

   
   plistdf %>%
    # make NA Stage if not the first sample
    ggplot(., aes(Pop1, Pop2, color = Stage)) + 
    geom_point(alpha = 0.9, colour= c(NA, "#333333")) +    # change Stage so colors don't show up
    theme_bw(base_size = 20) + 
    facet_grid(MigRate~Rep) + 
    theme(strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 90), # change the x-axis txt, too big
      #legend.position = "none", 
      legend.position="bottom", 
      legend.title = element_blank(), 
      panel.spacing = unit(0.1, "lines")) 

  # save plots
  ggsave("Rand_Sample_LogN_Normalised_Pop1vPop2.pdf", width = 210, height = 297, units = "mm")
} # end function 

# This makes the dataframe needed for P1 v P2 scatter plots 
plotNV = function(data_in, numReps){ # x is the migration rate
  out = as.data.frame(NULL)   
  for (i in 1:length(numReps)){
    r = numReps[i]
    data = data_in[[1]]
    xN= sapply(list(data[[r]][[1]]$TotalPopSize[c(START):c(END)], data[[r]][[2]]$TotalPopSize[START:END]), function(x) log10(x)) # get the log10 of N
    xN= as.data.frame(cbind(rep(r, nrow(xN)), seq(1:nrow(xN)), xN, c(rep("Before", Collect_Reps), rep("After", c(total_cycles-c(BurnIn_Reps+Collect_Reps))))))
    colnames(xN)<- c("Rep", "Time", "Pop1", "Pop2", "Stage")
    xNOut = xN %>% 
      mutate(., Stage = factor(Stage, levels = c("Before", "After"))) %>%
      mutate_at("Pop1",as.character) %>%
      mutate_at("Pop2",as.character) %>%
      mutate_at("Pop1",as.numeric) %>%
      mutate_at("Pop2",as.numeric) %>%
      tidyr::pivot_longer(Pop1:Pop2, names_to = "Pops", values_to= "LogN") %>%
      mutate(NormLogN = normalize(LogN)) #from Helper_Functions.R
    out = rbind(out, xNOut)  
  } # end 
  return(as.data.frame(out))
} #end function

# This gets population Fst, makes plots, and saves data for a sub-sample of population replicates at different migration rates
getFstData= function(input, N = N, numRepList = numRep_List, lambdaList){ 
  # Fst is only being measured here using the population allele frequency 
  x<- lapply(1:length(lambdaList), function(i) lapply(1:numReps, function(j) SamplePopsFS(input[[i]], numRep_List, j)))
  # combine the lists of lambda
  # collapse the list 
  plist = lapply(1:length(lambdaList), function(i){ 
    lapply(1:length(x[[1]]), function(j){
    y = as.data.frame(x[[i]][[j]])
    name = lambdaList[[i]]
    y %>%
      tibble::add_column(., MigRate = rep(name, nrow(y))) %>%
      mutate_at(vars(Loci1, Loci2, Loci3, Loci4, Loci5, Loci6, Loci7, Loci8, Loci9, Loci10), ~as.numeric(as.character(.))) %>%
      mutate_at(vars(Time), ~as.numeric(as.character(.))) %>%
      dplyr::mutate(Mean = rowMeans(.[,2:11]))
    })
  })
  #get all data into one structure
  plistdf = as.data.frame(do.call(rbind, unlist(plist, recursive = FALSE)))
  #write.table(plistdf, file= "Fst_Data.tsv", sep = "\t", row.names = F) # save data
  
  # get the Fst before migration 
  
  # find the time when Fst after migration == 50% of pre-migration Fst
  preFst =    plistdf %>%
    dplyr::filter(Time == c(BurnIn_Reps + Collect_Reps)) %>%
    dplyr::select(MigRate, Rep, Mean) %>%
    dplyr::rename(FstBefore = Mean)
  stablepostFst = plistdf %>%
    dplyr::filter(.,Time %in% c(BurnIn_Reps + Collect_Reps + 1:migrationEnd)) %>%
    dplyr::group_by(., MigRate, Rep) %>%
    dplyr::summarise(FstAfter = mean(Mean))
  # tibble with values where Fst after migration == 50% of pre-migration Fst
  fst50 = left_join(preFst, stablepostFst, by = c("MigRate", "Rep")) %>%
    dplyr::mutate(., halfFst = (FstBefore+FstAfter)/2) %>% # get the halfway point between premig fst and post mig equalibrium fst 
    
  #write.table(fst50, file = "Fst_50.tsv", sep = "\t", row.names = F)
  # make a boxplot
  fst50 %>%
    mutate_at(., "MigRate", as.factor) %>%
    ggplot(.,aes(x=MigRate, y=FstBefore, group=as.factor(MigRate))) + 
    geom_boxplot() + 
    theme_bw(base_size = 20)
  
ggsave("mean-premigFst.pdf", width = 210, height = 297, units = "mm")
  
  
# find the first time point where fst reaches 50% of pre and post fst
left_join(plistdf, fst50, by = c("MigRate", "Rep")) %>%
    dplyr::filter(.,Time %in% c(BurnIn_Reps + Collect_Reps + 1:migrationEnd)) 
    dplyr::filter(Mean <= halfFst) %>%
    dplyr::group_by(.,MigRate, Rep) %>%
    dplyr::filter(row_number() == 1) %>%
    ungroup() %>%
    mutate_at(., "MigRate", as.factor) %>%
    ggplot(.,aes(x=MigRate, y=Time, group=as.factor(MigRate))) + 
    geom_boxplot() + 
    ylim(c(BurnIn_Reps + Collect_Reps + 1), NA) + 
    labs(y = "post-Migration Fst = Fst50", x = "Migration Rate") +
    theme_bw(base_size = 20)
# 
ggsave("post-Migration-Fst50.pdf", width = 210, height = 297, units = "mm")

  
} # end function 
  

makeFstPlot= function(input, N = N, numRepList = numRep_List, lambdaList, subSample = subSample){ # subsample is the number of reps to sample
  # Fst is only being measured here using the population allele frequency 
  x<- lapply(1:length(lambdaList), function(i) lapply(subSample, function(j) SamplePopsFS(input[[i]], numRep_List, j)))
  # combine the lists of lambda
  # collapse the list 
  plist = lapply(1:length(lambdaList), function(i){ 
    lapply(1:length(x[[1]]), function(j){
      y = as.data.frame(x[[i]][[j]])
      name = lambdaList[[i]]
      y %>%
        tibble::add_column(., MigRate = rep(name, nrow(y))) %>%
        mutate_at(vars(Loci1, Loci2, Loci3, Loci4, Loci5, Loci6, Loci7, Loci8, Loci9, Loci10), ~as.numeric(as.character(.))) %>%
        mutate_at(vars(Time), ~as.numeric(as.character(.))) %>%
        dplyr::mutate(Mean = rowMeans(.[,2:11]))
    })
  })
  #get all data into one structure
  plistdf = as.data.frame(do.call(rbind, unlist(plist, recursive = FALSE)))
  plistdf %>%
    ggplot(., aes(Time, Mean, group = 1)) +
    geom_line(alpha = 0.9) +
    theme_bw(base_size = 20) +
    scale_x_continuous(breaks=seq(min(plistdf$Time)-1, max(plistdf$Time)-1, 100)) +
    facet_grid(MigRate~Rep) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_text(angle = 90), # change the x-axis txt, too big
          #legend.position = "none",
          legend.position="bottom",
          legend.title = element_blank(),
          panel.spacing = unit(0.1, "lines")) +
    scale_color_grey()
  # save plots
  ggsave("Rand_Sample_Fst_Pop1vPop2.pdf", width = 210, height = 297, units = "mm")
  } # end function plot fst

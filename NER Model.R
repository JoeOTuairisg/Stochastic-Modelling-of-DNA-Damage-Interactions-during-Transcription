library(magrittr)

projectRoot <- "/PATH TO SAVE OUTPUTS TO/"

saveplot = TRUE

set.seed(1)

geneLen <- 4*10^4 # gene length
drate <- 0.0025 # damage rate for starting damage amounts, damage/bp
drate <- 0.00049 # damage rate for starting damage amounts, damage/bp
meanDam <- geneLen*drate # expected number of damages on gene at start of sim
rnapSpeed <- 2000 # bp / min, how quickly RNAP transcribes
globalRepairRate <- 0.00034 # damages / min (repaired)
pd <- 0.9
pr <- 0.3 
fragDecay <- 0.0693 # as their half-lives is ~10 minutes, this is the rate from k = ln(2)/t_half

# original state of gene is ON
gene_state <- 1

#TSS Location
TSS <- 5000

# Number of RNAPs starting on gene immediately
nPocketRNAPs <- 1


# Params
ksyn <- 0.001639829 
off_on <- 0.006610024
on_off <- 0.0002912391
pd <- 0.9
pr <- 0.3
fragDecay <- 0.0693 
rnapSpeed <- 3500
drate <- 0.00049
globalRepairRate <- 0.00034

fnDistNextDamageRNAP <- function(damagePos,relRnapPos){
  return(damagePos[damagePos>relRnapPos][1]-relRnapPos) # first get a list of damage positions that come after where RNAP currently is
}

fnRNAPrates <- function(damagePos, relRnapPos, pd, pr, rnapSpeed){
  deltaDist <- fnDistNextDamageRNAP(damagePos,relRnapPos) #distance to next damage
  timeDeltaD <- deltaDist/rnapSpeed # how much time to get to next piece of damage
  rateDeltaD <- 1/timeDeltaD # the rate at which the next damage is encountered i.e. damage encountered / minute
  
  
  # getting the propensities for the event that occurs for the current damage, not probabilities yet
  # this will be used to find the time until the next event, regardless of what that event will be
  rateDetectRestart <- rateDeltaD*pd*pr # encounters damage, detects damage, restarts
  rateDetectNoRestart <- rateDeltaD*pd*(1-pr) # encounters damage, detects damage, doesn't restart
  rateNoDetect <- rateDeltaD*(1-pd) # encounters damage, doesn't detect damage
  
  return(c(rateDetectRestart, rateDetectNoRestart, rateNoDetect))
  
}

fnDoSimFixedTime = function(damagePos, rnapPos, repairFragPos, gene_state, 
                            tcur, tfin, ksyn,
                            pd, pr, fragDecay, rnapSpeed, drate_gen){
  
  # damagePos is locations on genes of damages
  # rnaPos is locations of where RNAP's are currently bound
  # repairFragPos is the locations at which repair has occurred, before they have been degraded
  # gene_state is 1 if gene is ON, 0 if gene is OFF
  # tcur is the time in the simulation that we are at
  # tfin is the last timepoint we simulate to
  # ksyn is the rate at which new RNAP's bind to the genes TSS
  # pd probability of damage detection for the RNAPs
  # pr probability of restarting transcription, having 
  # fragDecay is fragment decay rate
  # rnapSpeed is bp/min of RNAP moving along DNA
  
  
  
  while (tcur<tfin){
    
    if (length(damagePos) == 0){
      rnapPos <- numeric(0)
    } else {
      rnapPos <- rnapPos[rnapPos < max(damagePos)] # don't care about RNAP at last damage - can't produce any more frags
    }
    # for every RNAP, find the rates at which all of them 
    # encounters damage, detects damage, restarts
    # encounters damage, detects damage, doesn't restart
    # encounters damage, doesn't detect damage
    # for the closest damage in front of them
    perRNAPeventRates <- lapply(rnapPos, function(x) fnRNAPrates(damagePos,x,pd,pr,rnapSpeed)) # now each index of perRNAPeventRates is linked with the RNAPs in rnapPos
    
    # now we need to calculate the propensity of the rate that an event will happen to each RNAP    
    eachRNAPRate <- sapply(perRNAPeventRates, function(x) sum(x)) # eachRNAPRate is thus also linked with rnaPos via indices
    
    if (length(eachRNAPRate)>0 & gene_state == 1){ # checking if RNAPs are bound and the gene is in an ON state
      totalRnapRate <- sum(eachRNAPRate) # rate at which an RNAP event happens
    } else {
      totalRnapRate = 0 # no RNAPs on DNA, nothing will happen
    }
    
    nDamages = length(damagePos) # getting amount of damages on gene for calculating global repair
    if (nDamages == 0) {
      globalRepair = 0
    }
    else {
      globalRepair = nDamages * globalRepairRate
    }
    
    # setting the gene switching rate depending on if the gene is active or not
    gene_switch_rate <- ifelse(gene_state == 1, on_off, off_on)
    
    nRepFrag <- length(repairFragPos)
    
    if (gene_state == 0) {
      ksyn_curr <- 0
    } else {
      ksyn_curr <- ksyn
    }
    
    nextEventRates = c(ksyn_curr, #new RNAP initiates
                       nRepFrag*fragDecay, # rate of fragments broken down 
                       globalRepair, # global repair event
                       gene_switch_rate, # gene turns on / off
                       drate_gen * geneLen, # damage accumulates on gene
                       totalRnapRate) # something happens to one of the RNAPs
    
    deltaT = rexp(1,rate = sum(nextEventRates)) # sample exponential dist to find time until next event (unknown still what event)
    tcur = tcur+deltaT # progressing forward in time, for while loop
    
    if (tcur>tfin){ # past the end point of simulation
      break()
    }
    # perform 1 draw, once. We make a vector of probabilities (i.e. normalise them so they sum to 1, 
    # but are still proportional). e.g. next nexteventRates = [2, 0.16, 5.72], then normalise and
    # choose an index with probability = whatever is at that index
    # This returns a vector of the chosen event like (0, 1, 0) and by
    # using which == 1 here we find out which event was chosen, by selecting the index
    whichEvent <- which(rmultinom(1,1,nextEventRates/sum(nextEventRates))==1)
    
    # Logic for the different events
    if (whichEvent == 1){ # RNAP initiates
      start_site <- TSS
      rnapPos = c(start_site, rnapPos) # add a new RNAP to the very first position of the gene
      
    } else if (whichEvent == 2) { #repair fragment decays
      
      # Making a vector of probabilities, of length nRepFrag, with each of them having the same probability of being chosen
      whichFrag <- which(rmultinom(1,1,rep(1/nRepFrag,nRepFrag))==1)
      repairFragPos = repairFragPos[-whichFrag] # removing the repair fragment for that genic position
      nRepFrag = length(repairFragPos) # update the amount of excised fragments present
      
    } else if (whichEvent == 3) { # global repair event, entirely separate from RNAP activity
      whichDam <- which(rmultinom(1,1,rep(1/nDamages, nDamages))==1) # choosing a damage index randomly from the vector of damage positions
      globalDamagePos <- damagePos[whichDam] # finding out where that damage is on the genome
      
      # repairing and creating a repair fragment
      repairFragPos =  c(repairFragPos, globalDamagePos)
      # updating amount of repair fragments present in system
      nRepFrag = length(repairFragPos)
      
      # repair the damage, remove from damagePos vector
      damagePos = damagePos[-whichDam]
      
    } else if (whichEvent == 4) { # gene switches state event
      gene_state <- ifelse(gene_state == 1, 0, 1)
      
    } else if (whichEvent == 5) { # Damage Accumulates
      
      newDamage <- runif(1, 0, geneLen) # finding the position of the new damagse
      damagePos <- c(damagePos, newDamage) %>% sort() # add to damagePos vector and sort
      
      
    } else if (whichEvent == 6) { # RNAP event
      
      
      # we need to sample our eachRNAPRate to find out which of the RNAPs the event will occur to
      whichRNAP = which(rmultinom(1,1, eachRNAPRate/sum(eachRNAPRate))==1)
      
      # now we need to find out what event occurs to that RNAP, so get the vector of propensities for that RNAP
      relRates <- perRNAPeventRates[[whichRNAP]]
      # now sample the rates of the chosen RNAP to find what event occurs to this RNAP
      whichRNAPEvent = which(rmultinom(1,1, relRates/sum(relRates))==1)
      
      if (whichRNAPEvent == 1){ # rnap moves, detects, restarts
        beforePos = rnapPos[whichRNAP] # finding where the chosen RNAP started at, before moving 
        stallingDamagePos = damagePos[damagePos>beforePos][1] # finding the gene position of the damage we have encountered
        
        # create a repair fragment at that location, update amount of repair fragments in system
        repairFragPos =  c(repairFragPos,stallingDamagePos)
        nRepFrag = length(repairFragPos)
        
        # "move" the RNAP
        # remove RNAP from the positions vector
        rnapPos = rnapPos[-whichRNAP]
        # add RNAP back at stalling damage
        rnapPos = c(rnapPos, stallingDamagePos)
        
        # remove damage from damages vector
        damagePos = damagePos[!(damagePos %in% c(stallingDamagePos))]
        
        
      } else if (whichRNAPEvent == 2){ # rnap moves, detects, no restarts
        
        # position rnap was at before this event occurred (moving to site of damage)
        beforePos = rnapPos[whichRNAP]
        stallingDamagePos = damagePos[damagePos>beforePos][1] # position of next damage
        
        # repairing and creating a repair fragment 
        repairFragPos =  c(repairFragPos,stallingDamagePos)
        # updating amount of repair fragments present in system
        nRepFrag = length(repairFragPos)
        
        # RNAP dissociates, doesn't reattach to gene
        rnapPos = rnapPos[-whichRNAP]
        # but it did repair the damage, so remove this
        damagePos = damagePos[!(damagePos %in% c(stallingDamagePos))]
        
      } else if(whichRNAPEvent == 3) { # rnap moves, no detection
        
        # RNAP position before event
        beforePos = rnapPos[whichRNAP]
        stallingDamagePos = damagePos[damagePos>beforePos][1] # location of damage
        
        # RNAP moves, but no repair fragment is created as no detection occurred
        # so remove old position from rnapPos
        rnapPos = rnapPos[-whichRNAP]

        # add new position to rnapPos
        rnapPos = c(rnapPos, stallingDamagePos)
        
      }
      
      
      
    }
    
  }
  
  return(list(damagePos, rnapPos, repairFragPos, gene_state))
  
}

fnSimRepairFragPos <- function(tDetections, # give the timepoints we wish to sample at
                               ksyn,pd,pr, # then give system parameters
                               fragDecay,rnapSpeed,geneLen,
                               drate, drate_gen){
  
  
  deltaTDetect = diff(c(0,tDetections))
  
  simOutsTimePoints <- list()
  
  # initial number of damages on gene
  nInitDam <- rpois(1,geneLen*drate)
  # Uniform Damage Pattern
  damagePos <- runif(nInitDam, 0, geneLen) %>% sort()

    rnapPos <- c(TSS)
  repairFragPos = c()
  gene_state <- 1
  
  
  simOutsTimePoints[[1]] = list(damagePos, rnapPos, repairFragPos, gene_state)
  
  
  for ( k in 1:length(tDetections)){
    
    relSims = fnDoSimFixedTime(simOutsTimePoints[[k]][[1]], # damagePos
                               simOutsTimePoints[[k]][[2]], # rnapPos
                               simOutsTimePoints[[k]][[3]], # repairFragPos
                               simOutsTimePoints[[k]][[4]], # gene_state
                               0,
                               deltaTDetect[k],
                               ksyn,
                               pd,
                               pr,
                               fragDecay,
                               rnapSpeed, 
                               drate_gen)
    simOutsTimePoints[[k+1]] = relSims
  }
  
  
  # returning damage Map, RNAP position map, repair Map
  damageMap <- lapply(2:length(simOutsTimePoints), function(j) simOutsTimePoints[[j]][[1]])
  rnapMap <- lapply(2:length(simOutsTimePoints), function(j) simOutsTimePoints[[j]][[2]])
  repairMap <- lapply(2:length(simOutsTimePoints), function(j) simOutsTimePoints[[j]][[3]])
  
  
  maps <- list("damageMap" = damageMap, "rnapMap" = rnapMap, "repairMap" = repairMap)
  return(maps)
  
}


tcur = 0 # current time, starting at 0
nRepFrag = length(repairFragPos) # number of repair fragments

# this gives the outputs at each time point for each interval
simOutsTimePoints <- list()

# initialise rnaps on gene to 0
rnapPos <- c()

# initialise repair fragment count to 0
repairFragPos = c()

set.seed(1)
nsim = 10000
tDetections = c(1, 4, 8, 16, 24) * 60
deltaTDetect = diff(c(0,tDetections))
resultsRepairFragPos <- list()
resultsDamageFragPos <- list()

library(doParallel)
library(parallel)
numCores <- parallel::detectCores(logical = FALSE) - 1

cl <- makeCluster(numCores)
registerDoParallel(cl)


start_time <- proc.time()

# this is running the simulation nsim times, thus returning nsim repairFragPos 
# vectors at each time point in tDetections in parallel
# for example then this is returning 5000 lists of length 4
simOutMultPos = foreach(x = 1:nsim, .packages = c("magrittr")) %dopar% {
  # fnSimRepairFragPos(tDetections, ksyn, pd, pr, fragDecay, rnapSpeed, geneLen, drate)
  #fnSimRepairFragPos(tDetections,ksyn=.2,pd=.4,pr=.1, 1/100, .25*2000, 2*10^4, 2/(2*10^4))
  fnSimRepairFragPos(tDetections, ksyn, pd, pr, fragDecay, rnapSpeed, geneLen, drate, drate_gen)
}

stopCluster(cl)

# Getting the repair Map for every time point from every simulation
simOutMultRepairFragPosTimeStrat <- lapply(1:length(tDetections), function(k){
  sapply(1:nsim, function(j) simOutMultPos[[j]]$repairMap[[k]]) %>% unlist %>% 
    return()
})

# Getting the damage Map for every time point from every simulation
simOutMultDamageFragPosTimeStrat <- lapply(1:length(tDetections), function(k){
  sapply(1:nsim, function(j) simOutMultPos[[j]]$damageMap[[k]]) %>% unlist %>% 
    return()
})

# Getting the rnapPosition Map for every time point from every simulation
simOutMultRnapPosTimeStrat <- lapply(1:length(tDetections), function(k){
  sapply(1:nsim, function(j) simOutMultPos[[j]]$rnapMap[[k]]) %>% unlist %>% 
    return()
})

# for more detailed locations, split up xgrid into more points
xgrid = seq(0,geneLen,length.out = 41)

# here we are going over every repair fragment present at each of our time points
# we are summing up how many of them are in between certain locations on the gene, as determined by xgrid
fragCntsGrid = sapply(1:length(tDetections), function(k){
  sapply(1:(length(xgrid)-1), function(j){
    sum(simOutMultRepairFragPosTimeStrat[[k]]>=xgrid[j] & # summing over True and Falses is just 1s and 0s 
          simOutMultRepairFragPosTimeStrat[[k]]<xgrid[j+1]) %>% # so this line is counting the instances of each repair frag being less than a certain threshold on the gene
      return()
  }) %>% return() 
}) %>% t()

damCntsGrid = sapply(1:length(tDetections), function(k){
  sapply(1:(length(xgrid)-1), function(j){
    sum(simOutMultDamageFragPosTimeStrat[[k]]>=xgrid[j] & # summing over True and Falses is just 1s and 0s 
          simOutMultDamageFragPosTimeStrat[[k]]<xgrid[j+1]) %>% # so this line is counting the instances of each damage being less than a certain threshold on the gene
      return()
  }) %>% return() 
}) %>% t()

rnapPosGrid = sapply(1:length(tDetections), function(k){
  sapply(1:(length(xgrid)-1), function(j){
    sum(simOutMultRnapPosTimeStrat[[k]]>=xgrid[j] & # summing over True and Falses is just 1s and 0s 
          simOutMultRnapPosTimeStrat[[k]]<xgrid[j+1]) %>% # so this line is counting the instances of each rnap being less than a certain threshold on the gene
      return()
  }) %>% return() 
}) %>% t()

resultsRepairFragPos <- fragCntsGrid
resultsDamageFragPos <- damCntsGrid
resultsRnapPos <- rnapPosGrid

end_time <- proc.time()
print(end_time - start_time)

# midpoints of xgrid segments, just so we plot our counted number of repair fragments 
# in the middle of these gridpoints rather than the start, just for prettier plotting
xgridMids = (xgrid[-1]+xgrid[1:(length(xgrid)-1)])/2
cols = c("red","gold","green","cyan")


library(ggplot2)
library(tidyr)
library(dplyr)

if (saveplot == TRUE){ 
  
  repair_data_list <- list()
  damage_data_list <- list()
  
  
  # make RepairMap into a dataframe
  repair_results_df <- as.data.frame(resultsRepairFragPos)
  repair_results_df <- as.data.frame(t(repair_results_df))
  repair_results_df$Position <- xgridMids
  names(repair_results_df) <- c(paste0(tDetections, " min"), "Position")
  
  repair_pivoted_df <- pivot_longer(repair_results_df, paste0(tDetections, " min"), names_to="tDetect", values_to="RepairSites")
  
  repair_data_list <- repair_pivoted_df
  
  
  # make DamageMap into a dataframe
  damage_results_df <- as.data.frame(resultsDamageFragPos)
  damage_results_df <- as.data.frame(t(damage_results_df))
  damage_results_df$Position <- xgridMids
  names(damage_results_df) <- c(paste0(tDetections, " min"), "Position")
  
  damage_pivoted_df <- pivot_longer(damage_results_df, paste0(tDetections, " min"), names_to="tDetect", values_to="Damages")
  
  damage_data_list <- damage_pivoted_df
  
  # make rnapPos into a dataframe
  rnap_results_df <- as.data.frame(resultsRnapPos)
  rnap_results_df <- as.data.frame(t(rnap_results_df))
  rnap_results_df$Position <- xgridMids
  names(rnap_results_df) <- c(paste0(tDetections, " min"), "Position")
  
  rnap_pivoted_df <- pivot_longer(rnap_results_df, paste0(tDetections, " min"), names_to="tDetect", values_to="rnaps")
  
  rnap_data_list <- rnap_pivoted_df

    # reordering levels for plotting
  repair_data_list$tDetect <- factor(
    repair_data_list$tDetect,
    levels = paste0(tDetections, " min")
  )
  
  damage_data_list$tDetect <- factor(
    damage_data_list$tDetect,
    levels = paste0(tDetections, " min")
  )
  
  rnap_data_list$tDetect <- factor(
    rnap_data_list$tDetect,
    levels = paste0(tDetections, " min")
  )
  
  repair_plot <- ggplot(repair_data_list, aes(x = Position,
                                              y = RepairSites,
                                              color = factor(tDetect),
                                              group = tDetect)) +
    geom_line(linewidth=0.25) +
    geom_vline(xintercept = TSS, colour="darkgrey", linetype="dashed") +
    scale_x_continuous(
      name = "Position (bp)",
      breaks = c(pretty(repair_data_list$Position), TSS),
      labels = function(breaks) {
        sapply(breaks, function(x) if (x == TSS) "TSS" else as.character(x))
      }
    ) +
    labs(
      color = "Detection (min)",
      x     = "Position (bp)",
      y     = "Fragments detected"
    ) +
    theme_bw(base_size = 4) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.25, "cm"),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      plot.title = element_text(size = 21),
      axis.title = element_text(size=21),
      axis.text.x = element_text(size=16),
      axis.text.y = element_text(size=16)
    ) + 
    ggtitle(paste0("Repair Map"))
  
  
  filtered_rnap_df <- rnap_data_list %>% filter(Position >= 6000)
  
  ridgeline_plot <- ggplot(filtered_rnap_df, aes(x = Position, y = tDetect, fill=as.factor(tDetect))) + 
    geom_density_ridges(alpha=0.6, 
                        stat="binline", 
                        breaks=seq(0, geneLen, by = 1000), 
                        linewidth = 0,
                        color = "white") +
    coord_cartesian(xlim = c(0, 40000)) +
    scale_y_reverse(breaks = sort(unique(rnap_data_list$tDetect)),
                    labels = paste0(sort(unique(rnap_data_list$tDetect)), " min")) +
    theme_ridges() +
    theme(legend.position = "none",
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size=8)) +
    xlab("Gene Position") +
    ylab("Time of Detection")
  
  damage_plot <- ggplot(damage_data_list, aes(x = Position,
                                              y = Damages,
                                              color = factor(tDetect),
                                              group = tDetect)) +
    geom_line(linewidth=0.25) +
    geom_vline(xintercept = TSS) +
    geom_vline(xintercept = TSS, colour="darkgrey", linetype="dashed") +
    scale_x_continuous(
      name = "Position (bp)",
      breaks = c(pretty(repair_data_list$Position), TSS),
      labels = function(breaks) {
        sapply(breaks, function(x) if (x == TSS) "TSS" else as.character(x))
      }
    ) +
    labs(
      color = "Detection (min)",
      x     = "Position (bp)",
      y     = "Damages detected"
    ) +
    theme_bw(base_size = 4) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.25, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 12),
      axis.title = element_text(size=12),
      axis.text.x = element_text(size=10),
      axis.text.y = element_text(size=10)
    ) + 
    ggtitle(paste0("Damage Map"))
  
  library(patchwork)
  repair_damage_plot <- repair_plot / ridgeline_plot / damage_plot + plot_layout(ncol = 1, axis_titles="collect_x") & theme_light(base_size = 26)
  print(repair_plot)
  
  print("Time elapsed: ")
  print((end_time - start_time)["elapsed"])
  
  mmToInches <- function(mm){return(0.03937008*mm)}
  plot_filename <- paste0(projectRoot,"/repair_ridgeline_damage_plot.pdf")
  ggsave(plot_filename, width=mmToInches(200), height=mmToInches(200))
  
}


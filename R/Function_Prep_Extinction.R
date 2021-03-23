# This function prepares the data for the extinction function: 2 inputs abundance and biomass matrix 
# Arguments: abundance matrix: Taxa to site
#            biomass matrix: Taxa to site
#            replicate: summed, average or none (you must have only one record per taxa and per site)
#            correction: TRUE or FALSE, informs whether you want to complete missing abundance/biomass value with the average individual biomass
#            BPi: TRUE or FALSE, informs whether you want to import the Mi and Ri score for BPi calculation

prep.extinct <-  function(x, y, rep = "sum", correction = TRUE, BPi = TRUE) {

  ## package required
  require(tidyverse)
  
  ## preparing the data
  abn <- as.data.frame(x)
  colnames(abn)[1]<-"taxa"
  bio <- as.data.frame(y)
  colnames(bio)[1]<-"taxa"
  
  dimcheck <- (dim(abn) != dim(bio))
  if(any(dimcheck == T)){
    stop("Dimensions of abundance and biomass do not match, please check that taxa/site corresponds in both dataframes")
  }
  
  ###Data formatting***********************************************
  #****************************************************************
  if (rep == "sum"){
    # If we want the sum of rep per station
    # Dematrifying abundance & sum of taxa per site
    abndf<-abn %>%
      gather(site, A, 2:length(abn)) %>% 
      group_by(taxa, site) %>% 
      summarise (Ai = sum(A)) %>% 
      data.frame()
    
    # Ordering
    abndf<-abndf[with(abndf, order(taxa, site)), ]
    
    # Dematrifying biomass & sum of taxa per site
    biodf<-bio %>%
      gather(site, B, 2:length(bio)) %>% 
      group_by(taxa, site) %>% 
      summarise(Bi = sum(B)) %>% 
      data.frame()
    
    # Ordering
    biodf<-biodf[with(biodf, order(taxa, site)), ]
    
  } else if (rep == "mean"){
    # if we want the mean of rep per station
    # Dematrifying abundance & mean of taxa per site
    abndf<-abn %>%
      gather(site, A, 2:length(abn)) %>% 
      group_by(taxa, site) %>% 
      summarise (Ai = mean(A, na.rm=T)) %>% 
      data.frame()
    
    #Ordering
    abndf<-abndf[with(abndf, order(taxa, site)), ]
    
    #Dematrifying biomass & sum of taxa per site
    biodf<-bio %>%
      gather(site, B, 2:length(bio)) %>% 
      group_by(taxa, site) %>% 
      summarise(Bi = mean(B)) %>% 
      data.frame()
    
    #Ordering
    biodf<-biodf[with(biodf, order(taxa, site)), ]
    
  } else {
    
    # If there is already only one observation per station
    # Dematrifying abundance & mean of taxa per site
    abndf<-abn %>%
      gather(site, A, 2:length(abn)) %>% 
      rename(Ai = A) %>% 
      data.frame()
    
    #Ordering
    abndf<-abndf[with(abndf, order(taxa, site)), ]
    
    #Dematrifying biomass & sum of taxa per site
    biodf<-bio %>%
      gather(site, B, 2:length(bio)) %>% 
      rename(Bi = B) %>% 
      data.frame()
    
    #Ordering
    biodf<-biodf[with(biodf, order(taxa, site)), ]
  }
  
  #Combining abundance and biomass
  Bi<-biodf$Bi; ab_full<-cbind(abndf, Bi)
  
  test<-ab_full %>% 
    group_by(taxa, site) %>% 
    tally()
  
  if(any(test$n > 1)){
    stop("There are multiple records of the same taxon per site, please make sure there is only one record per taxon and per site by using the rep argument")
  }
 
  #Calculating the individual average biomass
  ab_tot<-ab_full %>% 
    group_by(taxa) %>% 
    summarise(Atot = sum(Ai),
              Btot = sum(Bi)) %>% 
    data.frame()
  ab_tot$Bind<-ab_tot$Btot/ab_tot$Atot
  
  #taxa-specific individual average biomass
  ab_full<-left_join(ab_full, ab_tot[, c("taxa", "Bind")]) 
  
  
  #Few missing data abundance value but no biomass & vice-versa
  abnmiss<-dim(ab_full[ab_full$Ai ==0 & ab_full$Bi != 0, ])[1]
  biomiss<-dim(ab_full[ab_full$Ai !=0 & ab_full$Bi == 0, ])[1]
  
  if (abnmiss >1 ){
    print(paste0("WARNING: There are ", abnmiss, " records with biomass > 0 and abundance == 0, if correction = T, estimates are generated with individual biomass"))
  }
  
  if (biomiss >1 ){
    print(paste0("WARNING: There are ", biomiss, " records with abundance > 0 and biomass == 0, if correction = T, estimates are generated with individual biomass"))
  }
  
  if (correction == TRUE){
    
    #Completing the dataset
    #Abundace missing -> Estimates based on average biomass
    ab_full[ab_full$Ai ==0 & 
              ab_full$Bi != 0,
            'Ai']<-ab_full[ab_full$Ai ==0 &
                             ab_full$Bi != 0,
                            'Bi']/ab_tot[(ab_full[ab_full$Ai ==0 &
                                                    ab_full$Bi != 0,
                                                       'taxa']),'Bind']
    
    #Biomass missing -> Estimates based on abundance and average biomass
    ab_full[ab_full$Ai !=0 &
              ab_full$Bi == 0 | is.na(ab_full$Bi),
            'Bi']<-ab_full[ab_full$Ai !=0 &
                                  ab_full$Bi == 0 | is.na(ab_full$Bi),
                                'Bind']*ab_full[ab_full$Ai !=0 &
                                                  ab_full$Bi == 0 | is.na(ab_full$Bi), 'Ai']
    
    
  }
  
  if (BPi == TRUE){
  #Formatting BPi parameters
  load("C:/Users/cg05/OneDrive - CEFAS/Science/Project Cefas - Other/Russia/Project/Task/Workshop/data/BPi.RData")
    
  BPi<-BPi %>% 
    select(Species_Updated_20112017, Mi_final, Ri_final) %>% 
    distinct(Species_Updated_20112017, .keep_all=T) %>% 
    rename(taxa = Species_Updated_20112017, Mi = Mi_final, Ri = Ri_final)
  
  #Combine BPi parameter
  ab_full<-left_join(ab_full, BPi)
  
  #Missing BPi
  Mimiss<-dim(ab_full[ab_full$Mi == 0 & is.na(ab_full$Mi) == T, ])[1]
  Rimiss<-dim(ab_full[ab_full$Ri == 0 & is.na(ab_full$Ri) == T, ])[1]
  
  if (Mimiss >1 ){
    print(paste0("WARNING: There are ", Mimiss, " records with Mi score missing, it will general NAs upon BPi calculation"))
  }
  
  if (Rimiss >1 ){
    print(paste0("WARNING: There are ", Rimiss, " records with Ri score missing, it will general NAs upon BPi calculation"))
  }
  
  ab_full$BPi<-sqrt(ab_full$Bi/ab_full$Ai)*ab_full$Mi * ab_full$Ri
  
  alldata<-ab_full %>% 
    filter (Ai > 0) %>% 
    group_by(taxa, Bind, Mi, Ri) %>% 
    summarise (Ai = mean(Ai, na.rm=T),
               Bi = mean(Bi, na.rm=T)) %>% 
    data.frame()
  } else {
    
    alldata<-ab_full %>% 
      filter (Ai > 0) %>% 
      group_by(taxa, Bind) %>% 
      summarise (Ai = mean(Ai, na.rm=T),
                 Bi = mean(Bi, na.rm=T)) %>% 
      data.frame()
  }
  alldata
}

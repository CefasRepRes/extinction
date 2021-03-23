# This function prepares the data for the extinction function: 2 inputs abundance and biomass matrix 
# Arguments: log.BPc: TRUE or FALSE, do you want the BPc scores to be plot on a log scale (default)

plot.extinct <-  function(x, log.BPc = T) {
  
  ## package required
  pkgs = c("tidyverse", "mgcv", "hexbin", "RColorBrewer")
  for(p in pkgs){
    if(!require(p, character.only = TRUE)) install.packages(p)
    library(p, character.only = TRUE)
  } 
  
  ## Colours
  rf <- colorRampPalette(rev(brewer.pal(7,'Spectral')))
  r <- rf(200)
  
  ## Data
  output<-as.data.frame(x)
  
  if (log.BPc == T) {
    gExtinct<-ggplot(output, aes(Nsp, log(BPc))) +
      geom_jitter(colour="grey",alpha=0.1)+
      geom_density_2d(colour='black', alpha=0.3, bins=20)+
      stat_density2d(aes(fill = ..density.., alpha=..density..), 
                     geom = "tile", contour = FALSE)+
      scale_fill_gradientn(colours=r)+
      theme_bw()+
      theme(legend.position="none")+
      stat_smooth(geom='line', alpha=0.5, se=TRUE, colour='red', linetype="twodash", size=1)
    
  } else {
    
    gExtinct<-ggplot(output, aes(Nsp, BPc)) +
      geom_jitter(colour="grey",alpha=0.1)+
      geom_density_2d(colour='black', alpha=0.3, bins=20)+
      stat_density2d(aes(fill = ..density.., alpha=..density..), 
                     geom = "tile", contour = FALSE)+
      scale_fill_gradientn(colours=r)+
      theme_bw()+
      theme(legend.position="none")+
      stat_smooth(geom='line', alpha=0.5, se=TRUE, colour='red', linetype="twodash", size=1)
  }
  print(gExtinct)
}

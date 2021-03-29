# This function performs an extinction sequence according to a predefined probability order
# dataframe must have the following column: taxa, Abundance (Ai), Biomass (Bi) and BPi score - Mi & Ri
# the function prep.extinct format the data the correct way based on 2 taxa-to-site abundance and biomass matrices
# Choose whether biomass compensation is active or not (TRUE / FALSE)
# comp <- FALSE (default), compensation is inactive; comp <- TRUE, compensation is active
# nsims sets the number of simulation (default is 500)
# extprob defines the extinction probability, default is random. Test changes.


#' Performs an extinction sequence
#'
#' performs an extinction sequence according to a predefined probability order
#'
#' @param comp compensation yes: comp = T; no: comp = F
#' @param nsims of simulation, nsims: numeric
#' @param extprob probability of extinction: extprob
#'
#' @return  List with how much biomass will be lost with the doomed species and Increase all biomass values by same proportion to maintain initial community biomass

#'
#' @author Clement Garcia, \email{broman@@wisc.edu}
#' @references \url{http://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#'
#' @examples
#' extinction(x, comp = F, nsims = 500, extprob = 1/nrow(x)
#'
#' @export
#'




extinction <-  function(x, comp = F, nsims = 500, extprob = 1/nrow(x)) {



  ## preparing the data
  alldata <- as.data.frame(x)
  ## Calculating BPc (According to the formula from Solan et al. 2004)
  alldata$BPi <- sqrt(alldata$Bi/alldata$Ai)*alldata$Mi * alldata$Ri
  ## Storing the number of species as a seperate variable
  nsp <- nrow(alldata)
  ## Calculating the total biomass at the start
  OrigTotalBiomass <- sum(alldata$Bi)
  # Add extra columns to record the abundances, biomass and extinction/compensation probabilities
  # during each simulation
  alldata$AiSim <- NA
  alldata$BiSim <- NA
  alldata$EPSim <- NA

  ## Set up output dataframe to record the results
  # sim           - the simulation number;
  # Nsp           - the number of species remaining in the community;
  # extsp         - the species that will go extinct at the next step;
  # bioext        - biomass that is lost dues to the species going extinct;
  # BPc           - the BPc for the remaining community;
  # Note that the extsp column contains the species that will go extinct
  # at the NEXT STEP, i.e. the first row records the BPc of the full
  # community and details who is going to go extinct next.

  output <- expand.grid(sim = 1:nsims, Nsp = nsp:1,
                        extsp=NA, bioext=NA, BPc=NA)

  # Sort the output dataframe for clarity
  output <- output[order(output$sim),]

  # With everything now set up we can run the simulations

  for (sim_count in 1:nsims){
    cat("sim_count: ", sim_count, "\n")

    # Reset abundances and probabilities for the current simulation
    alldata$AiSim <- alldata$Ai
    alldata$BiSim <- alldata$Bi
    alldata$EPSim <- extprob
    alldata$EPSim <- alldata$EPSim / sum(alldata$EPSim) # Normalise

    # Count down from full number of species to 1,
    # Removing species according to extinction probability

    for (sp_count in nsp:1)
    {
      # Calculate total community Bioturbation Potential
      BPc <- sum(alldata$AiSim * alldata$BPi, na.rm=T)

      # Store BPc results in the output dataframe
      output[output$sim == sim_count & output$Nsp== sp_count, "BPc"] <- BPc

      # Randomly pick a species to go extinct based on probability
      extinct <- which(cumsum(alldata$EPSim)>=runif(1))[1]

      # How much biomass will be lost with the doomed species
      biolost <- alldata[extinct,"BiSim"]

      # Call compensation function before we lose the extinct species
      if(comp==TRUE)
      {
        # Increase all biomass values by same proportion to maintain initial community biomass
        alldata$AiSim <- alldata$AiSim * (OrigTotalBiomass / (OrigTotalBiomass - biolost))

      }


      # Record ID of who has gone extinct
      output[output$sim == sim_count & output$Nsp == sp_count, "extsp"] <- as.character(alldata$taxa[extinct])
      # Record amount of compensation biomass
      output[output$sim == sim_count & output$Nsp == sp_count, "bioext"] <-biolost

      # Extinction happens! Set abundance/biomass and probability of doomed species to 0
      alldata[extinct,c("AiSim", "BiSim", 'EPSim')] <- 0

      # Extinction probability reset for newly compensating species
      alldata$EPSim <- alldata$EPSim / sum(alldata$EPSim)# Normalise
    }
  }
  output
}

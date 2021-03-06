# pod
R scripts for generating POD related figures.

# Requirements

The scripts variously require the following r packages: R2HTML for HTML output, sp for Spatial objects, and rgeos for Spatial functions (also grDevices for adjustcolor to make transparent areas (in core)), 

    install.packages("R2HTML")
    install.packages('sp') 
    install.packages('rgeos')

# Scripts

simplefigure.r produces a couple of simple diagrams, one comparison of atmospheric stability conditions with estimates for canine POD, another two different lateral range curves based on an exponential function.  

    source("simplefigure.r") 

sweepfigures.r produces a set of diagrams of parallel sweeps and random sweeps (conforming to the assumptions behind the exponential detection function of having different sweeps placed independently in the search area, each sweep being small relative to the total area searched, each sweep being long relative to lateral detection range, see Koopman, 1946 p.28) 
  
    source("sweepfigures.r")

detectioncurves.r plots various lateral range curves and runs simulations to graph detection functions for those lateral range cuves under paralell sweeps and random sweeps at different coverages.

    source("detectioncurves.r")

navigationerror.r runs a simple simulation of the effect of correlated navigation errors on detection functions for arbitrary lateral range curves.  

    source("navigationerror.r")


Koopman, B.O. 1946. Search and Screening. Operations Evaluation Group Report 56:1-172


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

detectioncurves.r plots various lateral range curves and runs simulations to graph detection functions for those lateral range cuves under paralell sweeps and random sweeps at different coverages.

    source("detectioncurves.r")

navigationerror.r runs a simple simulation of the effect of correlated navigation errors on detection functions for arbitrary lateral range curves.  

    source("navigationerror.r")


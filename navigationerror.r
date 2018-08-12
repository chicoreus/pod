require("R2HTML")# HTML output
require("grDevices") # for adjustcolor to make transparent areas.
require('sp')    # Spatial objects
require('rgeos') # Spatial functions

iterations<-10000
#iterations<-1000
#iterations<-100
nav_subjects <- 100
repeat_steps <- 100
repeat_steps <- 10

# Graph heighs and widths
squarePlotW<-1200
squarePlotH<-1200

plotW<-1500
plotH<-1200

fontScaling<-2.5

outputDirectory <- file.path("output_navigationerror")
if (!dir.exists(outputDirectory)){ dir.create(outputDirectory) }

segmentMinX<-0
segmentMaxX<-100
segmentMinY<-0
segmentMaxY<-100


naverr.env = new.env()
naverr.env$percent<-5
naverr.env$percentlength<-5

# Define some lateral range functions

expLRC<-function(x)   {exp(-abs(x)   ) }
exp10LRC<-function(x) {exp(-abs(x)^10) }
exp10NinteyLRC<-function(x) {exp(-abs(x)^10)*0.9 }
exp2LRC<-function(x) {exp(-abs(x)^2) }
defLRC<-function(x) { ifelse( x>=1 | x<=-1 ,0,1) }
# An inverse cube lateral range function
invcuLRC<-function(x) { 1-exp((-2*.1*10*(8*5))/(100*x^2))  }

# Lateral range curves that take a given effective lateral range (provide an ESW of approximately twice lr)
expLRCw<-function(x,lr)   {exp(-abs(x/lr)   ) }
exp10LRCw<-function(x,lr) {exp(-abs(x/lr)^10) }
exp10NinteyLRCw<-function(x,lr) {exp(-abs(x/lr)^10)*0.9 }
exp2LRCw<-function(x,lr) {exp(-abs(x/lr)^2) }
defLRCw<-function(x,lr) { ifelse( x>=lr | x<=-lr ,0,1) }
# An inverse cube lateral range function
invcuLRCw<-function(x,lr) { 1-exp((-2*.1*10*(8*5))/(100*(x/(lr/1.5))^2))  }

# function to obtain the probability of any one of a series of probabilities
pUnionAnyOneOf <- function(p){ 1 - prod(1-p) }

# Define some Sweep Models 

sweepRandom <- function(sweepcount, eswRange, segmentMinX=0, segmentMaxX=100, segmentMinY=0, segmentMaxY=100, doPlot=TRUE) { 
   sensorRuns <- list()
   searchLength <- 0
   for(iter in seq(1,sweepcount)) {
     segmentLength <- 15
     x1 <- runif(1,min=0,max=100)
     y1 <- runif(1,min=0,max=100)
     angle <- runif(1,min=0,max=2*pi)
     x2 <- x1 + (segmentLength * sin(angle))
     y2 <- y1 + (segmentLength * cos(angle))
     searchLength<-searchLength+segmentLength
     if (doPlot) { lines(c(x1,x2),c(y1,y2),col="blue") }
     slope <- - 1/((diff(c(y1,y2))/diff(c(x1,x2))))
     w <- eswRange
     #xes <- c(xm(x1,w,slope),xp(x1,w,slope),xp(x2,w,slope),xm(x2,w,slope))
     #yes <- c(ym(y1,w,slope),yp(y1,w,slope),yp(y2,w,slope),ym(y2,w,slope))
     #if (doPlot) { polygon(xes,yes,col=adjustcolor("gray",alpha.f=0.2))  }

     #search<-matrix(c(xes,yes),ncol=2,nrow=4)
     #sp <- Polygon(search,hole=as.logical(0))
     sensorRun <- Line(cbind(c(x1,x2),c(y1,y2)))
     sensorRunLines <- Lines(sensorRun,toString(iter))
     sensorRuns[iter]<-sensorRunLines
   }
   spSearch <- SpatialLines(sensorRuns,proj4string=CRS("+proj=utm +datum=WGS84"))
   spSweptArea <- gBuffer(spSearch,width=eswRange,byid=TRUE)
   spSweptAreaFlat <- gBuffer(spSearch,width=eswRange,byid=TRUE,capStyle="FLAT")  # Area only to right and left of sweep line.
   if (doPlot) {
       segment<-matrix(c(c(segmentMinX,segmentMinX,segmentMaxX,segmentMaxX),c(segmentMinY,segmentMaxY,segmentMaxY,segmentMinY)),ncol=2,nrow=4)
       sp <- Polygon(segment,hole=as.logical(0))
       spSegment <- SpatialPolygons(list(Polygons(list(sp),'Segment')),proj4string=CRS("+proj=utm +datum=WGS84"))
       spSearchedArea <- gIntersection(spSweptArea,spSegment,byid=TRUE)
       plot(spSearchedArea,col=adjustcolor("gray",alpha.f=0.2),add=TRUE)  
       spAreaNotSearched <- gDifference(spSegment,gUnionCascaded(spSearchedArea))
       plot(spAreaNotSearched,add=TRUE,col="lightskyblue")
   }
   return(pairlist(sweepLines=spSearch,eswArea=spSweptArea,eswAreaFlat=spSweptAreaFlat))
}

sweepParalell <- function(sweepcount, eswRange, segmentMinX=0, segmentMaxX=100, segmentMinY=0, segmentMaxY=100, doPlot=TRUE) {
   sensorRuns <- list()
   length <- segmentMaxX - segmentMinX
   sweepSpacing <- length/sweepcount
   searchLength <- 0
   for(iter in seq(1,sweepcount+1)) {
     segmentLength <- segmentMaxY - segmentMinY
     x1 <- (sweepSpacing * iter) - (sweepSpacing)
     y1 <- segmentMinY
     x2 <- (sweepSpacing * iter) - (sweepSpacing)
     y2 <- segmentMaxY
     searchLength<-searchLength+segmentLength
     if (doPlot) { lines(c(x1,x2),c(y1,y2),col="blue") }
     slope <- - 1/((diff(c(y1,y2))/diff(c(x1,x2))))
     w <- eswRange
     #xes <- c(xm(x1,w,slope),xp(x1,w,slope),xp(x2,w,slope),xm(x2,w,slope))
     #yes <- c(ym(y1,w,slope),yp(y1,w,slope),yp(y2,w,slope),ym(y2,w,slope))
     #if (doPlot) { polygon(xes,yes,col=adjustcolor("gray",alpha.f=0.2))  }

     #search<-matrix(c(xes,yes),ncol=2,nrow=4)
     #sp <- Polygon(search,hole=as.logical(0))
     sensorRun <- Line(cbind(c(x1,x2),c(y1,y2)))
     sensorRunLines <- Lines(sensorRun,toString(iter))
     sensorRuns[iter]<-sensorRunLines
   }
   spSearch <- SpatialLines(sensorRuns,proj4string=CRS("+proj=utm +datum=WGS84"))
   spSweptArea <- gBuffer(spSearch,width=eswRange,byid=TRUE)
   spSweptAreaFlat <- gBuffer(spSearch,width=eswRange,byid=TRUE,capStyle="FLAT")  # Area only to right and left of sweep line.
   if (doPlot) {
       segment<-matrix(c(c(segmentMinX,segmentMinX,segmentMaxX,segmentMaxX),c(segmentMinY,segmentMaxY,segmentMaxY,segmentMinY)),ncol=2,nrow=4)
       sp <- Polygon(segment,hole=as.logical(0))
       spSegment <- SpatialPolygons(list(Polygons(list(sp),'Segment')),proj4string=CRS("+proj=utm +datum=WGS84"))
       spSearchedArea <- gIntersection(spSweptArea,spSegment,byid=TRUE)
       plot(spSearchedArea,col=adjustcolor("gray",alpha.f=0.2),add=TRUE)
       spAreaNotSearched <- gDifference(spSegment,gUnionCascaded(spSearchedArea))
       if (!is.null(spAreaNotSearched)) { plot(spAreaNotSearched,add=TRUE,col="lightskyblue") }
   }
   return(pairlist(sweepLines=spSearch,eswArea=spSweptArea,eswAreaFlat=spSweptAreaFlat))
}

# Sweep model with paralell sweeps with a simple cumulative navigation error added
sweepParalellError <- function(sweepcount, eswRange, segmentMinX=0, segmentMaxX=100, segmentMinY=0, segmentMaxY=100, doPlot=TRUE) {
   sensorRuns <- list()
   length <- segmentMaxX - segmentMinX
   sweepSpacing <- length/sweepcount
   nextX  <- 0
   error <- naverr.env$percent
   sweepNo <- 0
   searchLength <- 0
   for(iter in seq(1,ceiling(sweepcount/2)+1)) {
     sweepNo <- sweepNo +1
     segmentLength <- segmentMaxY - segmentMinY
     x1 <- nextX
     if (x1<0) {  x1 <- -x1 } # don't allow for searches that entirely miss the search area.
     y1 <- segmentMinY
     stepErr <- runif(1,min=-error, max=error)
     nextX <- nextX + stepErr   # moving from x1 to x2 has navigation error.
     x2 <- nextX
     if (x2<0) {  x2 <- -x2 } # don't allow for searches that entirely miss the search area.
     y2 <- segmentMaxY
     searchLength<-searchLength+segmentLength
     if (doPlot) { lines(c(x1,x2),c(y1,y2),col="blue") }
     slope <- - 1/((diff(c(y1,y2))/diff(c(x1,x2))))
     w <- eswRange
     sensorRun <- Line(cbind(c(x1,x2),c(y1,y2)))
     sensorRunLines <- Lines(sensorRun,toString(sweepNo))
     sensorRuns[sweepNo]<-sensorRunLines
     sweepNo <- sweepNo +1
     segmentLength <- segmentMaxY - segmentMinY
     nextX <- nextX + sweepSpacing  # moving from x2 of one sweep to x2 of second sweep moves sweep spacing.
     x2 <- nextX
     y2 <- segmentMaxY
     stepErr <- runif(1,min=-error, max=error)
     nextX <- nextX + stepErr  # moving back from x2 to x1 has navigation error
     x1 <- nextX
     y1 <- segmentMinY
     nextX <- nextX + sweepSpacing  # moving from x1 of second sweep to x1 of next sweep moves sweep spacing.
     searchLength<-searchLength+segmentLength
     if (doPlot) { lines(c(x1,x2),c(y1,y2),col="blue") }
     slope <- - 1/((diff(c(y1,y2))/diff(c(x1,x2))))
     w <- eswRange
     sensorRun <- Line(cbind(c(x1,x2),c(y1,y2)))
     sensorRunLines <- Lines(sensorRun,toString(sweepNo))
     sensorRuns[sweepNo]<-sensorRunLines
   }
   spSearch <- SpatialLines(sensorRuns,proj4string=CRS("+proj=utm +datum=WGS84"))
   spSweptArea <- gBuffer(spSearch,width=eswRange,byid=TRUE)  # area with detection in all directions at the end of each sweep line.
   spSweptAreaFlat <- gBuffer(spSearch,width=eswRange,byid=TRUE,capStyle="FLAT")  # Area only to right and left of sweep line.
   if (doPlot) {
       segment<-matrix(c(c(segmentMinX,segmentMinX,segmentMaxX,segmentMaxX),c(segmentMinY,segmentMaxY,segmentMaxY,segmentMinY)),ncol=2,nrow=4)
       sp <- Polygon(segment,hole=as.logical(0))
       spSegment <- SpatialPolygons(list(Polygons(list(sp),'Segment')),proj4string=CRS("+proj=utm +datum=WGS84"))
       spSearchedArea <- gIntersection(spSweptArea,spSegment,byid=TRUE)
       plot(spSearchedArea,col=adjustcolor("gray",alpha.f=0.2),add=TRUE)
       spAreaNotSearched <- gDifference(spSegment,gUnionCascaded(spSearchedArea))
       if (!is.null(spAreaNotSearched)) { plot(spAreaNotSearched,add=TRUE,col="lightskyblue") }
   }
   return(pairlist(sweepLines=spSearch,eswArea=spSweptArea,eswAreaFlat=spSweptAreaFlat))
}

# Sweep model with paralell sweeps with a simple cumulative navigation error in both bearing and length added
sweepParalellError2 <- function(sweepcount, eswRange, segmentMinX=0, segmentMaxX=100, segmentMinY=0, segmentMaxY=100, doPlot=TRUE) {
   sensorRuns <- list()
   length <- segmentMaxX - segmentMinX
   sweepSpacing <- length/sweepcount
   nextX  <- 0
   error <- naverr.env$percent
   errorl <- naverr.env$percentlength
   errorLength <- errorl/length # error exppressed as a percent onf the segment length
   sweepNo <- 0
   for(iter in seq(1,ceiling(sweepcount/2)+1)) {
     randSweepLengthError <- runif(1, min=-errorLength*100, max=errorLength*100)
     sweepNo <- sweepNo +1
     segmentLength <- segmentMaxY - segmentMinY + randSweepLengthError
     x1 <- nextX
     if (x1<0) {  x1 <- -x1 } # don't allow for searches that entirely miss the search area.
     y1 <- segmentMinY
     stepErr <- runif(1,min=-error, max=error)
     nextX <- nextX + stepErr   # moving from x1 to x2 has navigation error.
     x2 <- nextX
     if (x2<0) {  x2 <- -x2 } # don't allow for searches that entirely miss the search area.
     y2 <- segmentMaxY + randSweepLengthError
     if (doPlot) { lines(c(x1,x2),c(y1,y2),col="blue") }
     slope <- - 1/((diff(c(y1,y2))/diff(c(x1,x2))))
     sensorRun <- Line(cbind(c(x1,x2),c(y1,y2)))
     sensorRunLines <- Lines(sensorRun,toString(sweepNo))
     sensorRuns[sweepNo]<-sensorRunLines
     sweepNo <- sweepNo +1
     segmentLength <- segmentMaxY - segmentMinY
     nextX <- nextX + sweepSpacing  # moving from x2 of one sweep to x2 of second sweep moves sweep spacing.
     x2 <- nextX
     y2 <- segmentMaxY + randSweepLengthError
     stepErr <- runif(1,min=-error, max=error)
     nextX <- nextX + stepErr  # moving back from x2 to x1 has navigation error
     x1 <- nextX
     y1 <- segmentMinY  # assumes y=0 is a marked baseline which is correctly returned to on each pair of sweeps.
     nextX <- nextX + sweepSpacing  # moving from x1 of second sweep to x1 of next sweep moves sweep spacing.
     if (doPlot) { lines(c(x1,x2),c(y1,y2),col="blue") }
     slope <- - 1/((diff(c(y1,y2))/diff(c(x1,x2))))
     sensorRun <- Line(cbind(c(x1,x2),c(y1,y2)))
     sensorRunLines <- Lines(sensorRun,toString(sweepNo))
     sensorRuns[sweepNo]<-sensorRunLines
   }
   spSearch <- SpatialLines(sensorRuns,proj4string=CRS("+proj=utm +datum=WGS84"))
   spSweptArea <- gBuffer(spSearch,width=eswRange,byid=TRUE)  # area with detection in all directions at the end of each sweep line.
   spSweptAreaFlat <- gBuffer(spSearch,width=eswRange,byid=TRUE,capStyle="FLAT")  # Area only to right and left of sweep line.
   if (doPlot) {
       segment<-matrix(c(c(segmentMinX,segmentMinX,segmentMaxX,segmentMaxX),c(segmentMinY,segmentMaxY,segmentMaxY,segmentMinY)),ncol=2,nrow=4)
       sp <- Polygon(segment,hole=as.logical(0))
       spSegment <- SpatialPolygons(list(Polygons(list(sp),'Segment')),proj4string=CRS("+proj=utm +datum=WGS84"))
       spSearchedArea <- gIntersection(spSweptArea,spSegment,byid=TRUE)
       plot(spSearchedArea,col=adjustcolor("gray",alpha.f=0.2),add=TRUE)
       spAreaNotSearched <- gDifference(spSegment,gUnionCascaded(spSearchedArea))
       if (!is.null(spAreaNotSearched)) { plot(spAreaNotSearched,add=TRUE,col="lightskyblue") }
   }
   return(pairlist(sweepLines=spSearch,eswArea=spSweptArea,eswAreaFlat=spSweptAreaFlat))
}

# Sweep model with paralell sweeps with a simple cumulative navigation error in both bearing and length with random wandering added
sweepParalellError3 <- function(sweepcount, eswRange, segmentMinX=0, segmentMaxX=100, segmentMinY=0, segmentMaxY=100, doPlot=TRUE) {
   sensorRuns <- list()
   length <- segmentMaxX - segmentMinX
   sweepSpacing <- length/sweepcount
   if (sweepcount==1) { sweepSpacing=sweepSpacing*.75 }
   nextX  <- 0
   error <- naverr.env$percent
   errorl <- naverr.env$percentlength
   errorLength <- errorl/length # error exppressed as a percent onf the segment length
   sweepNo <- 0
   for(iter in seq(1,ceiling(sweepcount/2)+1)) {
     randSweepLengthError <- runif(1, min=-errorLength*100, max=errorLength*100)
     sweepNo <- sweepNo +1
     segmentLength <- segmentMaxY - segmentMinY + randSweepLengthError
     x1 <- nextX
     if (x1<0) {  x1 <- -x1 } # don't allow for searches that entirely miss the search area.
     y1 <- segmentMinY
     stepErr <- runif(1,min=-error, max=error)
     nextX <- nextX + stepErr   # moving from x1 to x2 has navigation error.
     x2 <- nextX
     if (x2<0) {  x2 <- -x2 } # don't allow for searches that entirely miss the search area.
     y2 <- segmentMaxY + randSweepLengthError
     xsubstep <- -(x1 - x2)/6  # one sixth of the difference between x1 and x2, points would be on straight line between them
     midstepsx <- c(x1,runif(5,min=x1-error/3,max=x1+error/3),x2)  # random movement left/right of points along the sweep
     midstepsx <- midstepsx + c(0,xsubstep,xsubstep*2,xsubstep*3,xsubstep*4,xsubstep*5,0)  # so that deviation is off line from x1 to x2
     midstepsy <- c(y1,length/6,length/6*2,length/6*3,length/6*4,length/6*5,y2) # evenly spaced y positions of substeps
     midstepsy<-midstepsy + c(0,runif(5,min=-length/20,max=length/20),0) # add some random variation to y positions of substeps
     if (doPlot) { lines(midstepsx,midstepsy,col="blue") }
     sensorRun <- Line(cbind(midstepsx,midstepsy))
     sensorRunLines <- Lines(sensorRun,toString(sweepNo))
     sensorRuns[sweepNo]<-sensorRunLines
     sweepNo <- sweepNo +1
     segmentLength <- segmentMaxY - segmentMinY
     nextX <- nextX + sweepSpacing  # moving from x2 of one sweep to x2 of second sweep moves sweep spacing.
     x2 <- nextX
     y2 <- segmentMaxY + randSweepLengthError
     stepErr <- runif(1,min=-error, max=error)
     nextX <- nextX + stepErr  # moving back from x2 to x1 has navigation error
     x1 <- nextX
     y1 <- segmentMinY  # assumes y=0 is a marked baseline which is correctly returned to on each pair of sweeps.
     nextX <- nextX + sweepSpacing  # moving from x1 of second sweep to x1 of next sweep moves sweep spacing.
     midstepsx <- c(x1,runif(5,min=x1-error/3,max=x1+error/3),x2)
     midstepsy <- c(y1,length/6,length/6*2,length/6*3,length/6*4,length/6*5,y2)
     if (doPlot) { lines(midstepsx,midstepsy,col="blue") }
     sensorRun <- Line(cbind(midstepsx,midstepsy))
     sensorRunLines <- Lines(sensorRun,toString(sweepNo))
     sensorRuns[sweepNo]<-sensorRunLines
   }
   spSearch <- SpatialLines(sensorRuns,proj4string=CRS("+proj=utm +datum=WGS84"))
   spSweptArea <- gBuffer(spSearch,width=eswRange,byid=TRUE)  # area with detection in all directions at the end of each sweep line.
   spSweptAreaFlat <- gBuffer(spSearch,width=eswRange,byid=TRUE,capStyle="FLAT")  # Area only to right and left of sweep line.
   if (doPlot) {
       segment<-matrix(c(c(segmentMinX,segmentMinX,segmentMaxX,segmentMaxX),c(segmentMinY,segmentMaxY,segmentMaxY,segmentMinY)),ncol=2,nrow=4)
       sp <- Polygon(segment,hole=as.logical(0))
       spSegment <- SpatialPolygons(list(Polygons(list(sp),'Segment')),proj4string=CRS("+proj=utm +datum=WGS84"))
       spSearchedArea <- gIntersection(spSweptArea,spSegment,byid=TRUE)
       plot(spSearchedArea,col=adjustcolor("gray",alpha.f=0.2),add=TRUE)
       spAreaNotSearched <- gDifference(spSegment,gUnionCascaded(spSearchedArea))
       if (!is.null(spAreaNotSearched)) { plot(spAreaNotSearched,add=TRUE,col="lightskyblue") }
   }
   return(pairlist(sweepLines=spSearch,eswArea=spSweptArea,eswAreaFlat=spSweptAreaFlat))
}

# Function to calculate detection for some number of subjects in a segment under some lateral range curve and some sweep model.

calculateSweeps <- function(sweepcount,doPlot=TRUE, lrangecurve, sweepModel,  segmentMinX=0, segmentMaxX=100, segmentMinY=0, segmentMaxY=100,subjectCount=100) {

   if (doPlot) {  plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE) }
   searchLength <- 0
   # create a search segment 
   segment<-matrix(c(c(segmentMinX,segmentMinX,segmentMaxX,segmentMaxX),c(segmentMinY,segmentMaxY,segmentMaxY,segmentMinY)),ncol=2,nrow=4)
   sp <- Polygon(segment,hole=as.logical(0))
   spSegment <- SpatialPolygons(list(Polygons(list(sp),'Segment')),proj4string=CRS("+proj=utm +datum=WGS84"))

   searchEffort<- 0
   sweepWidth<-1

   # Determine the effective sweep range for the provided lateral range curve by integration 
   area<-integrate(function(x) { lrangecurve(x) },-100,100)
   eswRange = area$value/2

   # Create a search of the segment using the provided sweep model
   sweepRun <- sweepModel(sweepcount,eswRange,doPlot=doPlot)

   spSearch <- sweepRun$sweepLines
   # add in a half circle at the each end of each sweep to account for gDistance including distances off the ends of sweep lines, not just perpendiculars.
   spSweptArea <- sweepRun$eswArea
   spSweptAreaFlat <- sweepRun$eswAreaFlat  # with flat caps at the end of each sweep, so just perpendiculars.
   # Truncate the sweeps in the search at the bounaries of the search segment (search effort is limited to the search segment)
   spSearchedArea <- gIntersection(spSweptArea,spSegment,byid=TRUE,drop_lower_td=TRUE)
   spSearchSweeps <- gIntersection(spSearch,spSegment,byid=TRUE,drop_lower_td=TRUE)

   # Calculate the length of the path searched by the sweeps, and multipy this by the effective sweep width to obtain a measure of search effort as area searched.
   searchLength <- gLength(spSearchSweeps)
   searchEffort <- gArea(spSearchedArea)

   #spAreaSearched <- gIntersection(gUnionCascaded(spSearch),spSegment)  # area searched in segment, excluding search effort falling outside of segment

   # Place some (iterations) number of subjects in the search area
   subjectsx <- runif(subjectCount,min=segmentMinX,max=segmentMaxX)
   subjectsy <- runif(subjectCount,min=segmentMinY,max=segmentMaxY)
   spSubjects <- SpatialPoints(matrix(c(subjectsx,subjectsy),ncol=2),proj4string=CRS("+proj=utm +datum=WGS84"))
   subjects <- length(spSubjects)
 
   if(doPlot) { points(spSubjects,col="red",pch=19)    } 

   # calculate how many of the subjects are detected with the sweeps in the sweep model and the lateral range curve provided.
   detections <- 0
   for (x in 1:length(spSubjects)) { 
      # for each subject
      # determine the distances from each sweep to the subject
      distances<-gDistance(spSubjects[x],spSearchSweeps,byid=TRUE)

      detect<-FALSE
      detectP <- lrangecurve(distances)    # Find the detection probabilities for the subject at each distance
      toHit <- pUnionAnyOneOf(detectP)     # calculate the probability needed for detection of any one of the distances
      if (toHit >= runif(1)) { 
          # detect if a random number in the range 0-1 is less than the probability for any one detection.
          detect<-TRUE
          if(doPlot) { points(spSubjects[x],col="green",pch=19) } 
      }
      if (detect) { detections <- detections +1 } 
   }

   # Return the observed POD  as the number of detections per subject
   pod <- detections/subjects

   # spAreaNotSearched <- gDifference(spSegment,spAreaSearched)

   # Return coverage as the area searched over the area of the segment.
   coverage <- searchEffort / gArea(spSegment) 
   believedCoverage <- gArea(spSweptAreaFlat)/ gArea(spSegment)

   p<- 1 - exp(-( coverage ))  # for cross checking, calculate the POD for the exponential detection function with this coverage.
   #if (doPlot) { plot(spAreaNotSearched,add=TRUE,col="lightskyblue") }

   return( pairlist(p=p,coverage=coverage,podObs=pod,searchEffort=searchEffort,searchLength=searchLength,spSegment=spSegment,believedCoverage=believedCoverage) )
}

print("setup done")

# plots of model run

graphsweep<-"sweep200exp.png"
png(file.path(outputDirectory,graphsweep),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE) 

area<-integrate(function(x) { expLRC(x) },-100,100)
eswRange = area$value/2

#run <- calculateSweeps(200,doPlot=TRUE,lrangecurve=expLRC,sweepModel=sweepRandom,subjectCount=100)
sweepRandom(200,eswRange,doPlot=TRUE)

dev.off()

graphsweepsub<-"sweep200expsubjects.png"
png(file.path(outputDirectory,graphsweepsub),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE) 

area<-integrate(function(x) { expLRC(x) },-100,100)
eswRange = area$value/2

run <- calculateSweeps(200,doPlot=TRUE,lrangecurve=expLRC,sweepModel=sweepRandom,subjectCount=300)
#sweepRandom(200,eswRange,doPlot=TRUE)

dev.off()

# plots of model run - paralell sweeps

graphsweepp<-"sweepparalell30exp.png"
png(file.path(outputDirectory,graphsweepp),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE)

area<-integrate(function(x) { expLRC(x) },-100,100)
eswRange = area$value/2

#run <- calculateSweeps(200,doPlot=TRUE,lrangecurve=expLRC,sweepModel=sweepRandom,subjectCount=100)
sweepParalell(30,eswRange,doPlot=TRUE)

dev.off()

graphsweepsubp<-"sweepparalell30expsubjects.png"
png(file.path(outputDirectory,graphsweepsubp),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE)

area<-integrate(function(x) { expLRC(x) },-100,100)
eswRange = area$value/2

run <- calculateSweeps(30,doPlot=TRUE,lrangecurve=expLRC,sweepModel=sweepParalell,subjectCount=300)

dev.off()

naverr.env$percent<-5
naverr.env$percentlength<-5

graphsweepsubp<-"sweepparalell80naverrexp.png"
png(file.path(outputDirectory,graphsweepsubp),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE)

area<-integrate(function(x) { expLRC(x) },-100,100)
eswRange = area$value/2

sweepParalellError(80, eswRange, doPlot=TRUE)

dev.off()

graphsweepsubp2<-"sweepparalell80naverr2exp.png"
png(file.path(outputDirectory,graphsweepsubp2),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE)

area<-integrate(function(x) { expLRC(x) },-100,100)
eswRange = area$value/2

sweepParalellError2(80, eswRange, doPlot=TRUE)

dev.off()

graphsweepsubp3<-"sweepparalell80naverr3exp.png"
png(file.path(outputDirectory,graphsweepsubp3),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE)

fLRC<- function (x) { expLRCw(x,3) }
area<-integrate( fLRC,-100,100)
eswRange = area$value/2

sweepParalellError3(80, eswRange, doPlot=TRUE)

dev.off()



graphsweepsubp<-"sweepparalell80naverrexpsubj.png"
png(file.path(outputDirectory,graphsweepsubp),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE)

area<-integrate(function(x) { expLRC(x) },-100,100)
eswRange = area$value/2

run <- calculateSweeps(80,doPlot=TRUE,lrangecurve=expLRC,sweepModel=sweepParalellError,subjectCount=300)

dev.off()

graphsweepsubp<-"sweepparalell120naverrexp.png"
png(file.path(outputDirectory,graphsweepsubp),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE)

area<-integrate(function(x) { expLRC(x) },-100,100)
eswRange = area$value/2

sweepParalellError(120, eswRange, doPlot=TRUE)

dev.off()


graphsweepsubp<-"sweepparalell120naverrexpsubj.png"
png(file.path(outputDirectory,graphsweepsubp),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling)

plot(c(segmentMinX,segmentMaxX),c(segmentMinY,segmentMaxY),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE)

area<-integrate(function(x) { expLRC(x) },-100,100)
eswRange = area$value/2

run <- calculateSweeps(120,doPlot=TRUE,lrangecurve=expLRC,sweepModel=sweepParalellError,subjectCount=300)

dev.off()



print("random sweeps")

# setup for plotting
runsC<-c()
runsP<-c()

graphrs<-"randomsweepssim.png"
png(file.path(outputDirectory,graphrs),width=plotW,height=plotH)
par(cex=fontScaling)

# Sequence of numbers of sweeps to run to produce reasonably distributed points on detection function for the sweepRandom model
# over a coverage from about 0 to about 4
sweepnumbers <- c(seq(5,100,by=10),seq(120,200,by=20),seq(240,1000,by=40),seq(1040,1500,by=80))

for (sweepcount in sweepnumbers) { 
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=expLRC,sweepModel=sweepRandom,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsExpC<-runsC
runsExpP<-runsP
runsExpCRandom<-runsC
runsExpPRancom<-runsP
plot(runsC,runsP,xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection With Random Sweeps")
lines(runsC,runsP,lty=2)

runsC<-c()
runsP<-c()
for (sweepcount in sweepnumbers) { 
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=defLRC,sweepModel=sweepRandom,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsDefC<-runsC
runsDefP<-runsP
points(runsC,runsP,col="red")
lines(runsC,runsP,col="red",lty=2)

runsC<-c()
runsP<-c()
for (sweepcount in sweepnumbers) { 
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=exp10LRC,sweepModel=sweepRandom,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsExp10C<-runsC
runsExp10P<-runsP
points(runsC,runsP,col="darkgreen") 
lines(runsC,runsP,col="darkgreen",lty=2) 

runsC<-c()
runsP<-c()
for (sweepcount in sweepnumbers) { 
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=exp10NinteyLRC,sweepModel=sweepRandom,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsExp109C<-runsC
runsExp109P<-runsP
points(runsC,runsP,col="green") 
lines(runsC,runsP,col="green",lty=2) 

runsC<-c()
runsP<-c()
for (sweepcount in sweepnumbers) { 
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=invcuLRC,sweepModel=sweepRandom,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsIncC<-runsC
runsIncP<-runsP
points(runsC,runsP,col="green3") 
lines(runsC,runsP,col="green3",lty=2) 

# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
lines(runsC,f(runsC),col="blue")

legend("bottomright", pch = c(1, 1, 1, 1, 1, NA), lty=c(1,1,1,1,1,1), col = c("black", "red", "darkgreen","green", "green3", "blue"), 
        legend = c("Exponential LRC", "Definite Range LRC", "Exponential (10) LRC","exp(-abs(x)^10)*0.9 LRC", "Inverse Cube LRC", "Exponential Detection Function"))

dev.off()

# Plot just the exponential and definite range results
graphrsed<-"randomsweepssimexpdef.png"
png(file.path(outputDirectory,graphrsed),width=plotW,height=plotH)
par(cex=fontScaling)
plot(runsExpC,runsExpP,xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection With Random Sweeps")
lines(runsExpC,runsExpP,lty=2)
points(runsDefC,runsDefP,col="red")
lines(runsDefC,runsDefP,col="red",lty=2)

# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
lines(runsC,f(runsC),col="blue")

legend("bottomright", pch = c(1, 1, NA), lty=c(1,1,1), col = c("black", "red","blue"), 
        legend = c("Exponential LRC", "Definite Range LRC", "Exponential Detection Function"))

dev.off()


# Plot just the exponential and exponential ^10 * .9 results
graphrsee109<-"randomsweepssimexpexp109.png"
png(file.path(outputDirectory,graphrsee109),width=plotW,height=plotH)
par(cex=fontScaling)
plot(runsExpC,runsExpP,xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection With Random Sweeps")
lines(runsExpC,runsExpP,lty=2)
points(runsExp109C,runsExp109P,col="darkgreen")
lines(runsExp109C,runsExp109P,col="darkgreen",lty=2)

# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
lines(runsC,f(runsC),col="blue")

legend("bottomright", pch = c(1, 1, NA), lty=c(1,1,1), col = c("black", "darkgreen","blue"),
        legend = c("Exponential LRC", "exp(-abs(x)^10)*0.9 LRC", "Exponential Detection Function"))

dev.off()



rframe <- data.frame(runsExpC,runsExpP,runsDefC, runsDefP, runsExp10C, runsExp10P, runsIncC, runsIncP)
save(rframe,file=file.path(outputDirectory,paste('randomdetectionsim_',iterations,'.data',sep='')))

# Repeat but with paralell sweeps
print("paralell sweeps")

# Sequence of numbers of sweeps to run to produce reasonably distributed points on detection function for the sweepParalell model
# over a coverage from about 0 to about 4
sweepnumbers <- c(seq(1,20,by=2),seq(22,50,by=3),seq(52,200,by=5))

runsC<-c()
runsP<-c()
for (sweepcount in sweepnumbers) { 
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=expLRC,sweepModel=sweepParalell,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsExpC<-runsC
runsExpP<-runsP

runsC<-c()
runsP<-c()
for (sweepcount in sweepnumbers) { 
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=defLRC,sweepModel=sweepParalell,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsDefC<-runsC
runsDefP<-runsP

runsC<-c()
runsP<-c()
for (sweepcount in sweepnumbers) { 
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=exp10LRC,sweepModel=sweepParalell,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsExp10C<-runsC
runsExp10P<-runsP

runsC<-c()
runsP<-c()
for (sweepcount in sweepnumbers) {
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=exp10NinteyLRC,sweepModel=sweepParalell,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsExp109C<-runsC
runsExp109P<-runsP


runsC<-c()
runsP<-c()
for (sweepcount in sweepnumbers) { 
   run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=invcuLRC,sweepModel=sweepParalell,subjectCount=iterations)
   runsC <- c(runsC,run$coverage)
   runsP <- c(runsP,run$podObs)
}
runsIncC<-runsC
runsIncP<-runsP

# Plot detection functions for the exponential and definite range lateral range curves
graphpsed<-"paralellsweepssimexpdef.png"
png(file.path(outputDirectory,graphpsed),width=plotW,height=plotH)
par(cex=fontScaling)
plot(runsExpC,runsExpP,xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection With Paralell Sweeps")
lines(runsExpC,runsExpP,lty=2)
points(runsDefC,runsDefP,col="red")
lines(runsDefC,runsDefP,col="red",lty=2)
# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
lines(runsC,f(runsC),col="blue")
legend("bottomright", pch = c(1, 1, NA), lty=c(1,1,1), col = c("black", "red","blue"),
        legend = c("Exponential LRC", "Definite Range LRC", "Exponential Detection Function"))
dev.off()

# plot detection functions for exponential lateral range curve under random and parallel sweeps.

graphrpexpcomp<-"parallelrandomexpcomparison.png"
png(file.path(outputDirectory,graphpsed),width=plotW,height=plotH)
par(cex=fontScaling)
plot(runsExpC,runsExpP,xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Exponential Lateral Range Curve")
lines(runsExpC,runsExpP,lty=2)
points(runsExpCRandom,runsDefP,col="brown")
lines(runsExpCRandom,runsDefP,col="brown",lty=2)
legend("bottomright",pch=c(1,1,NA),lty=c(2,2,1),col=c("black","brown","blue"),legend=c("Parallel Sweeps","Random Sweeps","Exponential Detection Function"))
dev.off()

# Plot detection functions for the exponential, exponential with exponent of 10 scaled to 90%, and definite range lateral range curves
graphpsee109<-"paralellsweepssimexpe109.png"
png(file.path(outputDirectory,graphpsee109),width=plotW,height=plotH)
par(cex=fontScaling)
plot(runsExpC,runsExpP,xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection With Paralell Sweeps")
lines(runsExpC,runsExpP,lty=2)
points(runsDefC,runsDefP,col="red")
lines(runsDefC,runsDefP,col="red",lty=2)
points(runsExp109C,runsExp109P,col="darkgreen")
lines(runsExp109C,runsExp109P,col="darkgreen",lty=2)
# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
lines(runsC,f(runsC),col="blue")
legend("bottomright", pch = c(1,1,1, NA), lty=c(1,1,1,1), col = c("black","darkgreen","red","blue"),
        legend = c("Exponential LRC", "exp(-abs(x)^10)*0.9 LRC", "Definite Range LRC","Exponential Detection Function"))
dev.off()


# Plot detection functions for the exponential, exponential with exponent of 10, inverse cube, and definite range lateral range curves
graphps<-"paralellsweepssim.png"
png(file.path(outputDirectory,graphps),width=plotW,height=plotH)
par(cex=fontScaling)
plot(runsExpC,runsExpP,xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection With Paralell Sweeps")
lines(runsExpC,runsExpP,lty=2)
points(runsDefC,runsDefP,col="red")
lines(runsDefC,runsDefP,col="red",lty=2)
points(runsExp10C,runsExp10P,col="darkgreen") 
lines(runsExp10C,runsExp10P,col="darkgreen",lty=2) 
points(runsIncC,runsIncP,col="green3") 
lines(runsIncC,runsIncP,col="green3",lty=2) 
# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
lines(runsC,f(runsC),col="blue")
legend("bottomright", pch = c(1, 1, 1, 1, NA), lty=c(1,1,1,1,1), col = c("black", "red", "darkgreen", "green3", "blue"), 
        legend = c("Exponential LRC", "Definite Range LRC", "Exponential (10) LRC", "Inverse Cube LRC", "Exponential Detection Function"))
dev.off()

rframe <- data.frame(runsExpC,runsExpP,runsDefC, runsDefP, runsExp10C, runsExp10P, runsIncC, runsIncP,runsExp109C,runsExp109P)
save(rframe,file=file.path(outputDirectory,paste('paralelldetectionsim_',iterations,'.data',sep='')))


# Simulations with errors in navigation.  
# Model of navigation is of a sweep with some navigation error plus or minus of the intended end point, then a 
# fixed length step to the next sweep, then some navigation error plus or minus in returning to the base line, then 
# a fixed length step to the next sweep, then repeat.  Navigation errors can add up leaving gaps that aren't covered
# by other sweeps.  

print("parallel (with error) sweeps")

sweepnumbers <- c(seq(1,20,by=3),seq(22,50,by=4),seq(54,200,by=6))
sweepnumbers <- c(seq(1,10,by=2),seq(12,40,by=3),seq(43,150,by=6))

errorSteps <- seq(1,6)
runsExpCErr<-list()
runsExpBCErr<-list()
runsExpPErr<-list()
modelsExpErr<-list()
for (err in errorSteps) { 

   print(paste("Error: ",err))

   # Exponential lateral range curve with a lateral detection range of 3 (half the maximum error).
   fLRC<- function (x) { expLRCw(x,3) }
   naverr.env$percent<-err
   naverr.env$percentlength<-5
   runsC<-c()  # actual coverage
   runsBC<-c() # searchers believed coverage 
   runsP<-c()  # POD
   predictP<-0.8;
   for (sweepcount in sweepnumbers) { 
      rc <- list()
      rp <- list()
      mstep <- repeat_steps
      for (i in seq(1,mstep)) { 
          run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=fLRC,sweepModel=sweepParalellError,subjectCount=nav_subjects)
          rc <- c(rc,run$coverage)
          rp <- c(rp,run$podObs)
          rbc <- c(rc,run$believedCoverage)
      }
      rcp <- fr<-do.call(rbind, Map(data.frame, C=rc, P=rp,BC=rbc))
      rcpsorted <- rcp[with(rcp, order(C)),]
      runsC <- c(runsC, rcpsorted$C )  
      runsBC <- c(runsBC, rcpsorted$BC )  
      runsP <- c(runsP, rcpsorted$P )
   } 
   predictP<-mean(rcp$P)
   runsExp<- do.call(rbind, Map(data.frame, C=runsC, P=runsP,BC=runsBC))
   runsExp<-runsExp[order(runsExp$C),]
   runsExpCErr[[err]]<-runsExp$C
   runsExpBCErr[[err]]<-runsExp$BC
   runsExpPErr[[err]]<-runsExp$P
   m<-nls(P~a-exp(-abs(C+b)),data=runsExp,start=list(a=predictP,b=1))
   modelsExpErr[[err]]<-m
}

graphpcerr<-"paralellsweepserrsim.png"
png(file.path(outputDirectory,graphpcerr),width=plotW,height=plotH)
par(cex=fontScaling,lwd=2)

palette(gray(seq(.1,.6,len = length(errorSteps))))
palette(terrain.colors(length(errorSteps)))
palette(hsv(seq(.4,.1,len=length(errorSteps)),1,0.8,alpha=1))
plot(runsExpCErr[[1]],runsExpPErr[[1]],xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection (exponential LRC), Paralell Sweeps with navigation error")
for (err in errorSteps) { 
  points(runsExpCErr[[err]],runsExpPErr[[err]],col=palette()[err])
}
# make sure the lines overlay the points.
for (err in errorSteps) { 
  lines(runsExpCErr[[err]],predict(modelsExpErr[[err]]),col="black",lty=1)  
  lines(runsExpBCErr[[err]],predict(modelsExpErr[[err]]),col="firebrick4",lty=1) # searchers believed coverage, including navigating out of area.
  lines(runsExpCErr[[err]],predict(modelsExpErr[[err]]),col=palette()[err],lty=2)
}
# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
coverages<-seq(0,4,by=0.1)
lines(coverages,f(coverages),col="blue")

legends<-c()
pches<-c()
ltys<-c()
for (err in errorSteps) {  
   legends<-c(legends,paste("Nav Error%:",err," Goodn.Fit:",format(round( cor(runsExpCErr[[err]],predict(modelsExpErr[[err]])) ,2), nsmall=2) ))
   pches<-c(pches,1)
   ltys<-c(ltys,2)
}
legends<-c(legends,"Perceived vs Actual Coverage ")
pches<-c(pches,NA)
ltys<-c(ltys,1)
legends<-c(legends,"Exponential Detection Function")
pches<-c(pches,NA)
ltys<-c(ltys,1)
legend("bottomright", pch=pches, lty=ltys, col = c(palette(),"firebrick4","blue"), legend = legends )

palette("default")

dev.off()

# Plot again, but just with lines

graphpcerrf<-"paralellsweepserrsimfit.png"
png(file.path(outputDirectory,graphpcerrf),width=plotW,height=plotH)
par(cex=fontScaling,lwd=2)

palette(gray(seq(.1,.6,len = length(errorSteps))))
plot(runsExpCErr[[1]],runsExpPErr[[1]],col="white",xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection (exponential LRC), Paralell Sweeps with navigation error")
for (err in errorSteps) {
  #points(runsExpCErr[[err]],runsExpPErr[[err]],col=palette()[err])
  lines(runsExpCErr[[err]],predict(modelsExpErr[[err]]),col=palette()[err],lty=2)
}
# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
coverages<-seq(0,4,by=0.1)
lines(coverages,f(coverages),col="blue")

legends<-c()
pches<-c()
ltys<-c()
for (err in errorSteps) {
   m<-modelsExpErr[[err]]
   legends<-c(legends,paste("Nav Error%:",err," P~a-e^(-|C+b|) a:",format(round( coef(m)[[1]] ,2), nsmall=2 ),  " b:",format(round( coef(m)[[2]]  ,2), nsmall=2 ), "Goodn.Fit:",format(round( cor(runsExpCErr[[err]],predict(modelsExpErr[[err]])) ,2), nsmall=2) ) )
   pches<-c(pches,NA)
   ltys<-c(ltys,2)
}
legends<-c(legends,"Exponential Detection Function")
pches<-c(pches,NA)
ltys<-c(ltys,1)
legend("bottomright", pch=pches, lty=ltys, col = c(palette(),"blue"), legend = legends )

palette("default")

dev.off()

# again, but without the caption
graphpcerrfnc<-"paralellsweepserrsimfitnc.png"
png(file.path(outputDirectory,graphpcerrfnc),width=plotW,height=plotH)
par(cex=fontScaling,lwd=2.5)
palette(gray(seq(.1,.6,len = length(errorSteps))))
plot(runsExpCErr[[1]],runsExpPErr[[1]],col="white",xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection (exponential LRC), Paralell Sweeps with navigation error")
for (err in errorSteps) {
  #points(runsExpCErr[[err]],runsExpPErr[[err]],col=palette()[err])
  lines(runsExpCErr[[err]],predict(modelsExpErr[[err]]),col=palette()[err],lty=2)
}
# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
coverages<-seq(0,4,by=0.1)
lines(coverages,f(coverages),col="blue")
palette("default")

dev.off()


rframe <- data.frame(runsExpCErr)
save(rframe,file=file.path(outputDirectory,paste('paralellerr_cvals_',length(errorSteps),'_',nav_subjects,'.data',sep='')))
rframe <- data.frame(runsExpPErr)
save(rframe,file=file.path(outputDirectory,paste('paralellerr_pvals_',length(errorSteps),'_',nav_subjects,'.data',sep='')))



print("parallel (with length and bearing error) sweeps")

sweepnumbers <- c(seq(1,20,by=3),seq(22,50,by=4),seq(54,200,by=6))
sweepnumbers <- c(seq(1,10,by=2),seq(12,40,by=3),seq(43,150,by=6))
sweepnumbers <- c(seq(1,10,by=2),seq(12,40,by=3),seq(43,127,by=6))
sweepnumbers <- c(seq(1,10,by=2),seq(12,40,by=3),seq(43,115,by=6))

errorSteps <- seq(1,6)
runsExpCErr<-list()
runsExpBCErr<-list()
runsExpPErr<-list()
modelsExpErr<-list()
for (err in errorSteps) { 

   # Exponential lateral range curve with a lateral detection range of 3 (half the maximum error).
   fLRC<- function (x) { expLRCw(x,3) }
   print(paste("Error: ",err))

   naverr.env$percent<-err
   runsC<-c()  # actual coverage
   runsBC<-c() # searchers believed coverage 
   runsP<-c()  # POD
   predictP<-0.8;
   for (sweepcount in sweepnumbers) { 
      rc <- list()
      rp <- list()
      mstep <- repeat_steps
      for (i in seq(1,mstep)) { 
          run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=fLRC,sweepModel=sweepParalellError2,subjectCount=nav_subjects)
          rc <- c(rc,run$coverage)
          rp <- c(rp,run$podObs)
          rbc <- c(rc,run$believedCoverage)
      }
      rcp <- fr<-do.call(rbind, Map(data.frame, C=rc, P=rp,BC=rbc))
      rcpsorted <- rcp[with(rcp, order(C)),]
      runsC <- c(runsC, rcpsorted$C )  
      runsBC <- c(runsBC, rcpsorted$BC )  
      runsP <- c(runsP, rcpsorted$P )
   } 
   predictP<-mean(rcp$P)
   runsExp<- do.call(rbind, Map(data.frame, C=runsC, P=runsP,BC=runsBC))
   runsExp<-runsExp[order(runsExp$C),]
   runsExpCErr[[err]]<-runsExp$C
   runsExpBCErr[[err]]<-runsExp$BC
   runsExpPErr[[err]]<-runsExp$P
   m<-nls(P~a-exp(-abs(C+b)),data=runsExp,start=list(a=predictP,b=1))
   modelsExpErr[[err]]<-m
}

# Plot just with lines

graphpcerrfdb<-"paralellsweepserrsimfitdb.png"
png(file.path(outputDirectory,graphpcerrfdb),width=plotW,height=plotH)
par(cex=fontScaling,lwd=2)

palette(gray(seq(.1,.6,len = length(errorSteps))))
plot(runsExpCErr[[1]],runsExpPErr[[1]],col="white",xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection (exponential LRC), Paralell Sweeps with distance+bearing error")
for (err in errorSteps) {
  #points(runsExpCErr[[err]],runsExpPErr[[err]],col=palette()[err])
  lines(runsExpCErr[[err]],predict(modelsExpErr[[err]]),col=palette()[err],lty=2)
}
# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
coverages<-seq(0,4,by=0.1)
lines(coverages,f(coverages),col="blue")

legends<-c()
pches<-c()
ltys<-c()
for (err in errorSteps) {
   m<-modelsExpErr[[err]]
   legends<-c(legends,paste("Nav Error%:",err," P~a-e^(-|C+b|) a:",format(round( coef(m)[[1]] ,2), nsmall=2 ),  " b:",format(round( coef(m)[[2]]  ,2), nsmall=2 ), "Goodn.Fit:",format(round( cor(runsExpCErr[[err]],predict(modelsExpErr[[err]])) ,2), nsmall=2) ) )
   pches<-c(pches,NA)
   ltys<-c(ltys,2)
}
legends<-c(legends,"Exponential Detection Function")
pches<-c(pches,NA)
ltys<-c(ltys,1)
legend("bottomright", pch=pches, lty=ltys, col = c(palette(),"blue"), legend = legends )

palette("default")

dev.off()

rframe <- data.frame(runsExpCErr)
save(rframe,file=file.path(outputDirectory,paste('paralellerrdb_cvals_',length(errorSteps),'_',nav_subjects,'.data',sep='')))
rframe <- data.frame(runsExpPErr)
save(rframe,file=file.path(outputDirectory,paste('paralellerrdb_pvals_',length(errorSteps),'_',nav_subjects,'.data',sep='')))



print("parallel (with length and bearing error and wandering) sweeps")

sweepnumbers <- c(seq(1,20,by=3),seq(22,50,by=4),seq(54,200,by=6))
sweepnumbers <- c(seq(1,10,by=2),seq(12,40,by=3),seq(43,100,by=6))
sweepnumbers <- c(seq(1,10,by=2),seq(12,40,by=3),seq(43,85,by=6))
sweepnumbers <- c(seq(1,10,by=2),seq(12,40,by=3),seq(43,73,by=6))

errorSteps <- seq(1,6)
runsExpCErr<-list()
runsExpBCErr<-list()
runsExpPErr<-list()
modelsExpErr<-list()
for (err in errorSteps) { 

   # Exponential lateral range curve with a lateral detection range of 3 (half the maximum error).
   fLRC<- function (x) { expLRCw(x,3) }
   print(paste("Error: ",err))

   naverr.env$percent<-err
   runsC<-c()  # actual coverage
   runsBC<-c() # searchers believed coverage 
   runsP<-c()  # POD
   predictP<-0.8;
   for (sweepcount in sweepnumbers) { 
      rc <- list()
      rp <- list()
      mstep <- repeat_steps
      for (i in seq(1,mstep)) { 
          run <- calculateSweeps(sweepcount,doPlot=FALSE,lrangecurve=fLRC,sweepModel=sweepParalellError3,subjectCount=nav_subjects)
          rc <- c(rc,run$coverage)
          rp <- c(rp,run$podObs)
          rbc <- c(rc,run$believedCoverage)
      }
      rcp <- fr<-do.call(rbind, Map(data.frame, C=rc, P=rp,BC=rbc))
      rcpsorted <- rcp[with(rcp, order(C)),]
      runsC <- c(runsC, rcpsorted$C )  
      runsBC <- c(runsBC, rcpsorted$BC )  
      runsP <- c(runsP, rcpsorted$P )
   } 
   predictP<-mean(rcp$P)
   runsExp<- do.call(rbind, Map(data.frame, C=runsC, P=runsP,BC=runsBC))
   runsExp<-runsExp[order(runsExp$C),]
   runsExpCErr[[err]]<-runsExp$C
   runsExpBCErr[[err]]<-runsExp$BC
   runsExpPErr[[err]]<-runsExp$P
   m<-nls(P~a-exp(-abs(C+b)),data=runsExp,start=list(a=predictP,b=1))
   modelsExpErr[[err]]<-m
}

# Plot just with lines

graphpcerrfdbw<-"paralellsweepserrsimfitdbw.png"
png(file.path(outputDirectory,graphpcerrfdbw),width=plotW,height=plotH)
par(cex=fontScaling,lwd=2)

palette(gray(seq(.1,.6,len = length(errorSteps))))
plot(runsExpCErr[[1]],runsExpPErr[[1]],col="white",xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection (exponential LRC), Paralell Sweeps with wandering error")
for (err in errorSteps) {
  #points(runsExpCErr[[err]],runsExpPErr[[err]],col=palette()[err])
  lines(runsExpCErr[[err]],predict(modelsExpErr[[err]]),col=palette()[err],lty=2)
}
# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
coverages<-seq(0,4,by=0.1)
lines(coverages,f(coverages),col="blue")

legends<-c()
pches<-c()
ltys<-c()
for (err in errorSteps) {
   m<-modelsExpErr[[err]]
   legends<-c(legends,paste("Nav Error%:",err," P~a-e^(-|C+b|) a:",format(round( coef(m)[[1]] ,2), nsmall=2 ),  " b:",format(round( coef(m)[[2]]  ,2), nsmall=2 ), "Goodn.Fit:",format(round( cor(runsExpCErr[[err]],predict(modelsExpErr[[err]])) ,2), nsmall=2) ) )
   pches<-c(pches,NA)
   ltys<-c(ltys,2)
}
legends<-c(legends,"Exponential Detection Function")
pches<-c(pches,NA)
ltys<-c(ltys,1)
legend("bottomright", pch=pches, lty=ltys, col = c(palette(),"blue"), legend = legends )

palette("default")

dev.off()

graphpcerrdbw<-"paralellsweepserrsimdbw.png"
png(file.path(outputDirectory,graphpcerrdbw),width=plotW,height=plotH)
par(cex=fontScaling,lwd=2)

palette(gray(seq(.1,.6,len = length(errorSteps))))
palette(terrain.colors(length(errorSteps)))
palette(hsv(seq(.4,.1,len=length(errorSteps)),1,0.8,alpha=1))
plot(runsExpCErr[[1]],runsExpPErr[[1]],xlim=c(0,4),ylim=c(0,1),xlab="Coverage",ylab="POD",main="Detection (exponential LRC), Paralell Sweeps with wandering error")
for (err in errorSteps) {
  points(runsExpCErr[[err]],runsExpPErr[[err]],col=palette()[err])
}
# make sure the lines overlay the points.
for (err in errorSteps) {
  lines(runsExpCErr[[err]],predict(modelsExpErr[[err]]),col="black",lty=1)
  lines(runsExpBCErr[[err]],predict(modelsExpErr[[err]]),col="firebrick4",lty=1) # searchers believed coverage, including navigating out of area.
  lines(runsExpCErr[[err]],predict(modelsExpErr[[err]]),col=palette()[err],lty=2)
}
# Add the exponential detection function to the plot
f<-function(x) { 1- exp(-abs(x))}
coverages<-seq(0,4,by=0.1)
lines(coverages,f(coverages),col="blue")

legends<-c()
pches<-c()
ltys<-c()
for (err in errorSteps) {
   legends<-c(legends,paste("Nav Error%:",err," Goodn.Fit:",format(round( cor(runsExpCErr[[err]],predict(modelsExpErr[[err]])) ,2), nsmall=2) ))
   pches<-c(pches,1)
   ltys<-c(ltys,2)
}
legends<-c(legends,"Perceived vs Actual Coverage ")
pches<-c(pches,NA)
ltys<-c(ltys,1)
legends<-c(legends,"Exponential Detection Function")
pches<-c(pches,NA)
ltys<-c(ltys,1)
legend("bottomright", pch=pches, lty=ltys, col = c(palette(),"firebrick4","blue"), legend = legends )

palette("default")

dev.off()



rframe <- data.frame(runsExpCErr)
save(rframe,file=file.path(outputDirectory,paste('paralellerrdb_cvals_',length(errorSteps),'_',nav_subjects,'.data',sep='')))
rframe <- data.frame(runsExpPErr)
save(rframe,file=file.path(outputDirectory,paste('paralellerrdb_pvals_',length(errorSteps),'_',nav_subjects,'.data',sep='')))




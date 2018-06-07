require("R2HTML")# HTML output
require("grDevices") # for adjustcolor to make transparent areas.
require('sp')    # Spatial objects
require('rgeos') # Spatial functions

# ******** Configuration Parameters  *************

# Base font scaling parameter for graphs
fontScaling<-2

# Number of iterations for monte carlo simulations
iterations<-10000
#iterations<-100

# Graph heighs and widths
squarePlotW<-1200
squarePlotH<-1200

plotW<-1500
plotH<-1200

# Location to place output
outputDirectory <- file.path("output_detectioncurves")

# ******** Begin *************
if (!dir.exists(outputDirectory)){ dir.create(outputDirectory) }

output <- HTMLInitFile(outputDirectory,filename="pod_output")
HTML.title("Modeling POD");

HTML("This is output from the run of an R program that examines the relationship between lateral range curves, detection functions, coverage, and POD.")

print("setup done")
# Draw a diagram of several randomly oriented sweeps.

cf <- function(slope) { 1/sqrt(1+slope^2) } 
sf <- function(slope) { slope/sqrt(1+slope^2) }
xm <- function (x,d,slope) { x - d* cf(slope) } 
xp <- function (x,d,slope) { x + d* cf(slope) } 
ym <- function (y,d,slope) { y - d* sf(slope) } 
yp <- function (y,d,slope) { y + d* sf(slope) } 

graphrs<-"randomsweeps.png"
png(file.path(outputDirectory,graphrs),width=squarePlotW,height=squarePlotH)
par(mfrow=c(2,2))
par(cex.main=1.8, cex.lab=0.1)

for (maxiter in c(2,4,8,16)) { 

   plot(c(0,1),c(0,1),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE)
   searchLength <- 0
   for(iter in seq(1,maxiter)) { 
      side <- sample(1:4,1)
      side2 <- sample(1:4,1)
      while (side2==side) { side2 <- round(runif(1)*3)+1 }
      pos1 <- runif(1,min=-.1,max=1.1)
      while (pos1 <0 | pos1>1) {  pos1 <- runif(1,min=-.1,max=1.1) } 
      pos2 <- runif(1,min=-1.,max=1.1) 
      while (pos2 <0 | pos2>1) {  pos2 <- runif(1,min=-1.,max=1.1) }
      if (side==1) { 
        x1 <- pos1
        y1<- -0.1
      }
      if (side==3) { 
        x1 <- pos1
        y1<- 1.1
      }
      if (side==2) { 
        x1 <- 1.1
        y1<- pos1
      }
      if (side==4) { 
        x1 <- -0.1
        y1<- pos1
      }
    
      if (side2==1) {
        x2 <- pos2
        y2<- -0.1
      }
      if (side2==3) {
        x2 <- pos2
        y2<- 1.1
      }
      if (side2==2) {
        x2 <- 1.1
        y2<- pos2
      }
      if (side2==4) {
        x2 <- -0.1
        y2<- pos2
      }
      segmentLength <- sqrt( (abs(x1-x2)^2) * ( abs(y1-y2)^2))
      searchLength<-searchLength+segmentLength
      lines(c(x1,x2),c(y1,y2),col="blue")
      slope <- - 1/((diff(c(y1,y2))/diff(c(x1,x2))))
      w <- 0.1
      xes <- c(xm(x1,w,slope),xp(x1,w,slope),xp(x2,w,slope),xm(x2,w,slope))
      yes <- c(ym(y1,w,slope),yp(y1,w,slope),yp(y2,w,slope),ym(y2,w,slope))
      polygon(xes,yes,col=adjustcolor("gray",alpha.f=0.2))
   }
   # multipy by 10 (detectionRange = 1 instead of 0.1
   searchWidth<-2
   searchArea<-10*10
   searchLength<-searchLength*10
   w<-(searchWidth*searchLength)/searchArea
   p<-1-exp(-w)
   w1<-1-((searchWidth*searchLength)/(maxiter*searchArea))
   p1<-1-(w1^maxiter)
   title(main=paste("p=",round(p,digits=2)))
}
dev.off();

print("random sweeps diagram done")
# draw some random sweeps which fit the expectations of the exponential detection function (placed randomly, short relative to search area, longer than effective sweep width, placed independently of each other)

plotRandomSweeps <- function(maxiter,doPlot=TRUE) { 

   if (doPlot) {  plot(c(0,100),c(0,100),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE) }
   searchLength <- 0
   segment<-matrix(c(c(0,0,100,100),c(0,100,100,0)),ncol=2,nrow=4)
   sp <- Polygon(segment,hole=as.logical(0))
   spSegment <- SpatialPolygons(list(Polygons(list(sp),'Segment')),proj4string=CRS("+proj=utm +datum=WGS84"))
   sps<-vector(mode="list",length=maxiter)
   searchEffort<- 0
   for(iter in seq(1,maxiter)) {
     segmentLength <- 15
     x1 <- runif(1,min=0,max=100)
     y1 <- runif(1,min=0,max=100)
     angle <- runif(1,min=0,max=2*pi)
     x2 <- x1 + (segmentLength * sin(angle))
     y2 <- y1 + (segmentLength * cos(angle))
   
     searchLength<-searchLength+segmentLength
     if (doPlot) { lines(c(x1,x2),c(y1,y2),col="blue") }
     slope <- - 1/((diff(c(y1,y2))/diff(c(x1,x2))))
     w <- 2
     xes <- c(xm(x1,w,slope),xp(x1,w,slope),xp(x2,w,slope),xm(x2,w,slope))
     yes <- c(ym(y1,w,slope),yp(y1,w,slope),yp(y2,w,slope),ym(y2,w,slope))
     if (doPlot) { polygon(xes,yes,col=adjustcolor("gray",alpha.f=0.2))  }
     search<-matrix(c(xes,yes),ncol=2,nrow=4)
     sp <- Polygon(search,hole=as.logical(0))
     sptemp<-list(sp)
     sp<- gIntersection(spSegment,SpatialPolygons(list(Polygons(sptemp,'Segment')),proj4string=CRS("+proj=utm +datum=WGS84")))
     searchEffort<- searchEffort+gArea(sp)
     if (iter==1) {
       spSearch <- sp
     } else {
        spSearch <- gUnion(spSearch,sp)
     }
     sps[iter]<-sp
   }
   spAreaSearched <- gIntersection(gUnionCascaded(spSearch),spSegment)  # area searched in segment, excluding search effort falling outside of segment

   subjectsx <- runif(iterations,min=0,max=100)
   subjectsy <- runif(iterations,min=0,max=100)
   spSubjects <- SpatialPoints(matrix(c(subjectsx,subjectsy),ncol=2),proj4string=CRS("+proj=utm +datum=WGS84"))   
   detections <- length(gIntersection(spAreaSearched,spSubjects))
   pod <- detections/iterations

   spAreaNotSearched <- gDifference(spSegment,spAreaSearched)
   coverage <- searchEffort / gArea(spSegment)
   searchWidth<-2
   searchArea<-100*100
   w<-(searchWidth*searchLength)/searchArea
   p<- 1 - exp(-( coverage ))
   w1<-1-((searchWidth*searchLength)/(maxiter*searchArea))
   p1<-1-(w1^maxiter)
   if (doPlot) { title(main=paste("coverage=",round(coverage,digits=2)," pod=",round(p,digits=2))) }
   if (doPlot) { plot(spAreaNotSearched,add=TRUE,col="lightskyblue") }

   return( pairlist(p=p,coverage=coverage,podObs=pod) )
}

plotParallelSweeps <- function(totalSweeps,doPlot=TRUE) { 

   if (doPlot) {  plot(c(0,100),c(0,100),xaxs = 'i',yaxs = 'i', xaxt='n', yaxt='n', ann=FALSE) }
   searchLength <- 0
   segment<-matrix(c(c(0,0,100,100),c(0,100,100,0)),ncol=2,nrow=4)
   sp <- Polygon(segment,hole=as.logical(0))
   spSegment <- SpatialPolygons(list(Polygons(list(sp),'Segment')),proj4string=CRS("+proj=utm +datum=WGS84"))
   sps<-vector(mode="list",length=maxiter)
   searchEffort<- 0
   searchWidth<-2
   sweepSpacing = 100/totalSweeps 
   for(sweep in seq(1,totalSweeps)) {
     
     x1 <- (sweep * sweepSpacing) - (searchWidth) - (sweepSpacing/2)
     y1 <- 0
     x2 <- (sweep * sweepSpacing) + (searchWidth) - (sweepSpacing/2)
     y2 <- 100
   
     searchLength<-searchLength+100
     if (doPlot) { lines(c(x1+searchWidth,x2-searchWidth),c(y1,y2),col="blue") }
     slope <- - 1
     w <- 2
     xes <- c(x1,x1,x2,x2)
     yes <- c(y1,y2,y2,y1)
     if (doPlot) { polygon(xes,yes,col=adjustcolor("gray",alpha.f=0.2))  }
     search<-matrix(c(xes,yes),ncol=2,nrow=4)
     sp <- Polygon(search,hole=as.logical(0))
     sptemp<-list(sp)
     sp<- gIntersection(spSegment,SpatialPolygons(list(Polygons(sptemp,'Segment')),proj4string=CRS("+proj=utm +datum=WGS84")))
     searchEffort<- searchEffort+gArea(sp)
     if (sweep==1) {
       spSearch <- sp
     } else {
        spSearch <- gUnion(spSearch,sp)
     }
     sps[iter]<-sp
   }
   spAreaSearched <- gIntersection(gUnionCascaded(spSearch),spSegment)  # area searched in segment, excluding search effort falling outside of segment

   subjectsx <- runif(iterations,min=0,max=100)
   subjectsy <- runif(iterations,min=0,max=100)
   spSubjects <- SpatialPoints(matrix(c(subjectsx,subjectsy),ncol=2),proj4string=CRS("+proj=utm +datum=WGS84"))   
   detections <- length(gIntersection(spAreaSearched,spSubjects))
   pod <- detections/iterations

   spAreaNotSearched <- gDifference(spSegment,spAreaSearched)
   coverage <- searchEffort / gArea(spSegment)
   searchArea<-100*100
   w<-(searchWidth*searchLength)/searchArea
   p<- 1 - exp(-( coverage ))
   w1<-1-((searchWidth*searchLength)/(maxiter*searchArea))
   p1<-1-(w1^maxiter)
   if (doPlot) { title(main=paste("coverage=",round(coverage,digits=2)," pod=",round(p,digits=2))) }
   if (doPlot && !is.null(spAreaNotSearched)) { plot(spAreaNotSearched,add=TRUE,col="lightskyblue") }

   return( pairlist(p=p,coverage=coverage,podObs=pod) )
}

graphrsf<-"randomsweepfrags.png"
png(file.path(outputDirectory,graphrsf),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)

plotRandomSweeps(50,TRUE)

dev.off();

graphrsf2<-"randomsweepfrags2.png"
png(file.path(outputDirectory,graphrsf2),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)

plotRandomSweeps(190,TRUE)

dev.off();

graphrsf3<-"randomsweepfrags3.png"
png(file.path(outputDirectory,graphrsf3),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)

plotRandomSweeps(500,TRUE)

dev.off();

graphps<-"parallelsweep.png"
png(file.path(outputDirectory,graphps),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)

plotParallelSweeps(8,TRUE)

dev.off();


graphrsrsdf<-"randomsweepdetectfunct.png"
png(file.path(outputDirectory,graphrsrsdf),width=plotW,height=plotH)
par(cex=fontScaling)

byStep<-25

searchBits <- seq(from=25,to=500,by=byStep)
pr<-c()
co<-c()
for (makeSearches in searchBits) {
   iterate<-5
   for (rep in 1:iterate) { 
      result <- plotRandomSweeps(makeSearches,FALSE)
      pr<-c(pr,result$podObs)
      co<-c(co,result$coverage)
   }
}

rgframe <- data.frame(pr,co)
save(rgframe,file=file.path(outputDirectory,paste('randomsweeps_',byStep,'_',iterations,'.data',sep='')))

plot(co,pr, xlim=c(0,3),ylim=c(0,1),main="Randomly placed sweeps and the EDF",xlab="Coverage", ylab="POD")

searchBits <- seq(from=3,to=50,by=1)
prp<-c()
cop<-c()
for (makeSearches in searchBits) {
   iterate<-3
   for (rep in 1:iterate) {
      result <- plotParallelSweeps(makeSearches,FALSE)
      prp<-c(prp,result$podObs)
      cop<-c(cop,result$coverage)
   }
}

points(cop,prp)

rgframep <- data.frame(prp,cop)
save(rgframep,file=file.path(outputDirectory,paste('parallelsweeps_',iterations,'.data',sep='')))

coverage <- seq(0,4,by=0.1)
# exponential - any lateral range curve, random sweeps.
f<-function(x) {exp(-abs(x))}
pod <- c()
for (i in coverage) {
  area<-integrate(f, 0,i)
  pod<-c(pod,area$value)
}
lines(coverage,pod,col="blue")

f<-function(x) { ifelse( x>=1 | x<=-1 ,0,1) }
podp <- c()
for (i in coverage) {
  area<-integrate(f, 0,i)
  podp<-c(podp,area$value)
}
lines(coverage,podp,col="blue")

dev.off();

graphrsrsdf1<-"randomsweepdetectfunct1.png"
png(file.path(outputDirectory,graphrsrsdf1),width=plotW,height=plotH)
par(cex=fontScaling)

plot(co,pr, xlim=c(0,3),ylim=c(0,1),main="Randomly placed perfect broom sweeps fall on the EDF",xlab="Coverage", ylab="POD")
points(cop,prp)
lines(coverage,pod,col="blue")
lines(coverage,podp,col="blue")
lines(c(0,1,1,1),c(.63,.63,.63,0),col="brown")

dev.off()

print("random plots done")

HTML("<p>We'll start with lateral range curves.  Lateral range curves represent the probability that a target at some distance from a sensor will be detected by the sensor.  The x axis of a lateral range curve is the distance from the sensor to a target.  The y axis of a lateral range curve is the instantaneous probability of detection of a target at that range.</p>")

#  Distributions
graph1<-"distributioncurves.png"
png(file.path(outputDirectory,graph1),width=plotW,height=plotH)
par(mfrow=c(2,2))
# Exponential
x<-seq(-3,3,0.01)
y<-exp(-abs(x))
plot(x,y,type="l",main="exponential e^(-|x|)",mfg=c(1,1))
area<-integrate(function(x) {exp(-abs(x))}, -100,100)
lines(c(area$value/2,area$value/2),c(0,1),col="blue",mfg=c(1,1))
lines(c(-area$value/2,-area$value/2),c(0,1),col="blue")

# Normal
y<-exp(-abs(x)^2)
plot(x,y,type="l",main="e^(-|x|^2)",mfg=c(1,2))
area<-integrate(function(x) {exp(-abs(x)^2)}, -100,100)
lines(c(area$value/2,area$value/2),c(0,1),col="blue",mfg=c(1,2))
lines(c(-area$value/2,-area$value/2),c(0,1),col="blue")

y<-exp(-abs(x)^0.5)
plot(x,y,type="l",main="e^(-|x|^0.5)",mfg=c(2,1))
area<-integrate(function(x) {exp(-abs(x)^0.5)}, -100,100)
lines(c(area$value/2,area$value/2),c(0,1),col="blue",mfg=c(2,1))
lines(c(-area$value/2,-area$value/2),c(0,1),col="blue")

# Approximating perfect broom
y<-exp(-abs(x)^10)
plot(x,y,type="l",main="e^(-|x|^10)",mfg=c(2,2))
area<-integrate(function(x) {exp(-abs(x)^10)}, -100,100)
lines(c(area$value/2,area$value/2),c(0,1),col="blue",mfg=c(2,2))
lines(c(-area$value/2,-area$value/2),c(0,1),col="blue")
dev.off();



# perfect broom
graphbroom<-"broom.png"
png(file.path(outputDirectory,graphbroom),width=squarePlotW,height=squarePlotH)
par(cex=.9,cex.lab=fontScaling,cex.axis=fontScaling,mar=c(6, 6, 2, 2)+0.1,lwd=2)
par(cex=fontScaling);
f<-function(x) { ifelse( x>=1 | x<=-1 ,0,1) }
y<-f(x)
plot(x,y,type="h",main="Perfect Broom (definite range)",xlab="Distance",ylab="POD",col="gray95")
lines(x,y,col="black")
area<-integrate(f, -100,100)
lines(c(area$value/2,area$value/2),c(0,1),col="blue",mfg=c(2,2))
lines(c(-area$value/2,-area$value/2),c(0,1),col="blue")
arrows(area$value/2,0.6,-area$value/2,0.6,col="blue",code=3)
text(c(0),c(0.6),pos=3,col="blue",labels="Sweep Width")
dev.off();

HTMLInsertGraph(graphbroom,file=output,Caption="Perfect broom (definite range) lateral range curve, showing sweep width.")

HTML("<p>One simple lateral range curve is a definite range curve. Think of this as a perfect broom sweeping a floor.  Every target within the definite range (the width of the broom) is detected.  Every target outside of that definite range is not detected.   The perfect broom has an evident sweep width - the width of the broom.</p>")

# exponential lateral range curve (not to be confused with the exponential detection function)
graphedf<-"exp_lrc.png"
png(file.path(outputDirectory,graphedf),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling);
f<-function(x) {exp(-abs(x))}
y<-f(x)
plot(x,y,type="h",main="Exponential Lateral Range Curve",xlab="Distance",ylab="POD",col="gray95")
lines(x,y,col="black")
area<-integrate(f, -100,100)
lines(c(area$value/2,area$value/2),c(0,1),col="blue",mfg=c(2,2))
lines(c(-area$value/2,-area$value/2),c(0,1),col="blue")
arrows(area$value/2,0.6,-area$value/2,0.6,col="blue",code=3)
text(c(0),c(0.6),pos=3,col="blue",labels="ESW")
dev.off();

HTMLInsertGraph(graphedf,file=output,Caption="Exponential Lateral Range Curve, showing effective sweep width")

HTML("<p>Another lateral range curve models detection as falling off exponentially with distance from the sensor.  An exponential function is p(x)= e^|x|.  The instantaneous probability of detection at some distance x is a function of the mathematical constant e to the (absolute value of) x power.  Targets near to the sensor are very likely to be detected, but detection falls of exponentially with distance.  This is a reasonable model of visual detection under many circumstances.   Do not confuse the exponential lateral range curve with the exponential detection function - the two are mathematically similar but model very different aspects of search.</p>")

HTML("<p>The exponential lateral range curve, unlike the perfect broom (definite range curve), doesn't have a clear edge.  Because detection drops off continously from the sensor, it lacks an obvious range for talking about how wide a sweep is.  There is, however, an easy mathematical way to model a sweep width for the exponential lateral range curve.  If we draw the area under the curve as a rectangle with a height of 1, the width of the rectangle also happens to be the distance at which the number of objects detected inside that distance equals the number of objects missed outside that distance.   This distance is termed the Effective Sweep Width (and is plotted in blue).</p>")

# inverse cube
graphedic<-"inversecube.png"
png(file.path(outputDirectory,graphedic),width=squarePlotW,height=squarePlotH)
par(cex=fontScaling);
# values for k,h,A,a,b which produce a reasonable inverse cube function in range x -3 to 3
# f<-function(x) { 1-exp((-2*k*h*A)/(v*x^2))  }
f<-function(x) { 1-exp((-2*.1*10*(8*5))/(100*x^2))  }
y<-f(x)
plot(x,y,type="h",main="Inverse Cube Lateral Range Curve",xlab="Distance",ylab="POD",col="gray95")
lines(x,y,col="black")
area<-integrate(f, -100,100)
lines(c(area$value/2,area$value/2),c(0,1),col="blue",mfg=c(2,2))
lines(c(-area$value/2,-area$value/2),c(0,1),col="blue")
arrows(area$value/2,0.6,-area$value/2,0.6,col="blue",code=3)
text(c(0),c(0.6),pos=3,col="blue",labels="ESW")
dev.off();

HTMLInsertGraph(graphedic,file=output,Caption="Inverse cube lateral range curve, showing effective sweep width.")

HTML("<p>The inverse cube lateral range curve is a standard model for detection at sea.  It reflects a model of a sensor looking down on the ocean surface and a target wake which falls into a rectangle on the ocean surface.  The sensor's ability to detect the wake depends on the height of the sensor, the distance from the sensor to the target's wake, the area of the rectangle enclosing the target wake, and factors such as sea surface state.  This particular curve is for a set of parameters that illustrate the rapidly dropping off shape of the curve in the same arbitrary distance scale used in the other lateral range curve diagrams above.  As with the exponential lateral range curve, the area under the curve is the same as the area of a rectangle with height of 1 and a width equal to the effective sweep width.   Integrating to obtain the area under a lateral range curve gives the effective sweep width (for a rectangle with a height of p=1, conveniently easy to calculate as we are working with probabilities from 0 to 1, so the height of the rectangle is 1).</p>")


HTMLInsertGraph(graph1,file=output,Caption="Comparison of 4 functions that can model lateral range curves.")

HTML("<p>The 4 curves above are all variations on the exponential lateral range curve with different parameters.  Each also shows an effective sweep width calculated by integrating to find the area under the curve.</p>")

# pod/coverage
graphpc<-"coveragepod.png"
png(file.path(outputDirectory,graphpc),width=plotW,height=plotH)
par(cex=.9,cex.lab=fontScaling,cex.axis=fontScaling,mar=c(6, 6, 2, 2)+0.1,lwd=2)
#perfect broom/Definite range lateral range curve, parallel sweeps.
f<-function(x) { ifelse( x>=1 | x<=-1 ,0,1) }
coverage <- seq(0,4,by=0.1)
pod <- c()
for (i in coverage) { 
  area<-integrate(f, 0,i)
  pod<-c(pod,area$value)
}
plot(coverage,pod,type="l",xlab="Coverage",ylab="POD")
coverage_pb<-coverage
pod_pb<-pod
# exponential - any lateral range curve, random sweeps.
f<-function(x) {exp(-abs(x))}
pod <- c()
for (i in coverage) { 
  area<-integrate(f, 0,i)
  pod<-c(pod,area$value)
}
lines(coverage,pod,col="blue")
coverage_exp<-coverage
pod_exp<-pod

dev.off();


graphpcs1<-"coveragecurves.png"
png(file.path(outputDirectory,graphpcs1),width=plotW,height=plotH)
par(cex=fontScaling);

f<-function(x,s) { exp(-abs(x-s))}
lp<-function(f,x,spacing) { 
    r<-c()
    offsets<-seq(-1,12)
    for (offset in offsets) { 
       r<-c(r,f(x,offset*spacing))
    }
    return(r)
}

x<-c(0,10)
y<-c(0,1)
plot(x,y,type="p",xlim=c(0,10),ylim=c(0,1),main="Exponential lateral range curve, parallel sweeps",col="white",xlab="Distance",ylab="POD")

# simple case, probability of any one of a list of probabilities.
pUnionAnyOneOf <- function(p){ 1 - prod(1-p) } 

# General implementations, probability of n of a list of probabilities: 
# See: https://stackoverflow.com/questions/35119185/probability-of-the-union-of-three-or-more-sets
cp <- function(p) 
{
  # slow
  ev <- do.call(expand.grid,replicate(length(p),0:1,simplify=FALSE))
  pe <- apply(ev,1,function(x) prod(p*(x==1)+(1-p)*(x==0)))
  tapply(pe,rowSums(ev),sum)
}
cp.quadratic <- function(p) {
  P <- matrix(0, nrow=length(p), ncol=length(p))
  P[1,] <- rev(cumsum(rev(p * prod(1-p) / (1-p))))
  for (i in seq(2, length(p))) {
    P[i,] <- c(rev(cumsum(rev(head(p, -1) / (1-head(p, -1)) * tail(P[i-1,], -1)))), 0)
  }
  c(prod(1-p), P[,1])
}

# Given an x value, a function which returns a y value as
# a probability, and a set of offsets by which the function 
# is shifted on the x axis, return the probability represented
# by the union of any one of the y probabilities.
funion <- function(x,f,offsets) { 
   yvals<-c()
   for (offset in offsets) { 
       yval<-f(x,offset)
       yvals<-c(yvals,yval)
   }
   y<-pUnionAnyOneOf(yvals)
   return(y)
} 


# Given the environment f.env containing a function to evaluate which
# takes two parameters, an x and an offset (a function which can be invoked
# by funion) f.env$f, and a vector of offsets f.env$offsets), and given
# a vector of integers, evaluate funion(each x, the function, and the offsets.
# return a single value.
# This function allows integration of offsets of a function with funion 
# (which can't be directly integrated as it takes too many parameters 
# for integrate()).
fToIntegrate<-function(xes) { 
  retval <- c()
  for (x in xes) { 
     retval <- c(retval,funion(x,f.env$f,f.env$offsets))
  }
  return(retval)
} 

coverages<-c()
pods<-c()
sweepsPerWidth<-c(1,2,3,4,5,6,7,8,9,10)
width<-10
f<-function(x) {exp(-abs(x))}
area<-integrate(f,-100,100)
sweepWidth<-area$value
f<-function(x,s) { exp(-abs(x-s))}
for (sp in sweepsPerWidth) { 
   spacing = width / sp
   offsets<-seq(-1,15)
   offsets<-offsets * spacing
   f.env = new.env()
   f.env$offsets<-offsets
   f.env$f<-f
   if (sp==5) {
      xit<-seq(0,width,by=0.1)
      lines(xit,fToIntegrate(xit),col="blue")
      for (offset in offsets) { 
         lines(xit,f(xit,offset))
      }
      polygon(c(spacing-(sweepWidth/2),spacing-(sweepWidth/2),spacing+(sweepWidth/2),spacing+(sweepWidth/2)),c(0,1,1,0),col=adjustcolor("wheat",alpha.f=0.3),border="wheat3")
   }
   pod<-integrate(fToIntegrate,lower=0,upper=width,subdivisions=1000L)
   pod<-pod$value/width

   coverage = (sweepWidth*sp)/(width)
   coverages<-c(coverages,coverage)
   pods<-c(pods,pod)
} 
dev.off()

HTMLInsertGraph(graphpcs1,file=output,Caption="Exponential lateral range curve, showing multiple parallel sweeps.  One effective sweep width shown in tan.  Cumulative probability in blue.  Each separate sweep in black.")

coveragesFirst <- coverages
podsFirst <- pods

print("exponential lrc, parallel sweeps done")
#graphpcs1a<-"coveragecurves1a.png"
#png(file.path(outputDirectory,graphpcs1a),width=1500,height=1200)
#par(cex=fontScaling);

# random placement of sweeps, iterate multiple times and average exponential detection 
coverages<-c()
pods<-c()
sweepsPerWidth<-c(1,2,3,4,5,6,7,8,9,10)
width<-10
f<-function(x) {exp(-abs(x))}
area<-integrate(f,-100,100)
sweepWidth<-area$value
f<-function(x,s) { exp(-abs(x-s))}
# for y=e^(-|x]), random sweeps should have pod/coverage relationship of pod= integral from 0 to coverage of e^(-|coverage|)
for (sp in sweepsPerWidth) {
  piter <-0
  citer <- 0
  iter <- iterations
  for (i in 1:iter) { 
     spacing = width / sp
     offsets<- runif(sp,0,width)
  
     f.env = new.env()
     f.env$offsets<-offsets
     f.env$f<-f
     pod<-integrate(fToIntegrate,lower=0,upper=width,subdivisions=1000L)
     pod<-pod$value/width
     coverage = (sweepWidth*sp)/(width)
     piter <- piter + pod
     citer <- citer + coverage
  }
  coverages<-c(coverages,citer/iter)
  pods<-c(pods,piter/iter)
}

coveragesIter <- coverages
podsIter <- pods

# random placement of sweeps, iterate multiple times and average perfect broom
coverages<-c()
pods<-c()
sweepsPerWidth<-c(1,2,3,4,5,6,7,8,9,10)
width<-10
f<-function(x) { ifelse( x>=1 | x<=-1 ,0,1) }
area<-integrate(f,-100,100)
sweepWidth<-area$value
f<-function(x,s) { exp(-abs(x-s))}
# for y=e^(-|x]), random sweeps should have pod/coverage relationship of pod= integral from 0 to coverage of e^(-|coverage|)
for (sp in sweepsPerWidth) {
  piter <-0
  citer <- 0
  iter <- iterations
  for (i in 1:iter) { 
     spacing = width / sp
     offsets<- runif(sp,0,width)
  
     f.env = new.env()
     f.env$offsets<-offsets
     f.env$f<-f
     pod<-integrate(fToIntegrate,lower=0,upper=width,subdivisions=1000L)
     pod<-pod$value/width
     coverage = (sweepWidth*sp)/(width)
     piter <- piter + pod
     citer <- citer + coverage
  }
  coverages<-c(coverages,citer/iter)
  pods<-c(pods,piter/iter)
}

coveragesIterPB <- coverages
podsIterPB <- pods


# Normal lateral range curve, parallel sweeps
graphpcs2<-"coveragecurves2.png"
png(file.path(outputDirectory,graphpcs2),width=plotW,height=plotH)
par(cex=fontScaling);
x<-c(0,10)
y<-c(0,1)
plot(x,y,type="p",xlim=c(0,10),ylim=c(0,1),main="Normal exp(-|x^2|)",col="white",xlab="Distance",ylab="POD")

coverages<-c()
pods<-c()
sweepsPerWidth<-c(1,2,3,4,5,6,7,8,9,10)
width<-10
f<-function(x) {exp(-abs(x^2))}
area<-integrate(f,-100,100)
sweepWidth<-area$value
f<-function(x,s) { exp(-abs((x-s)^2)) }
for (sp in sweepsPerWidth) {
   spacing = width / sp
   offsets<-seq(-1,15)
   offsets<-offsets * spacing
   f.env = new.env()
   f.env$offsets<-offsets
   f.env$f<-f
   if (sp==5) {
      xit<-seq(0,width,by=0.1)
      lines(xit,fToIntegrate(xit),col="blue")
      for (offset in offsets) {
         lines(xit,f(xit,offset))
      }
      polygon(c(spacing-(sweepWidth/2),spacing-(sweepWidth/2),spacing+(sweepWidth/2),spacing+(sweepWidth/2)),c(0,1,1,0),col=adjustcolor("wheat",alpha.f=0.3),border="wheat3")
   }
   pod<-integrate(fToIntegrate,lower=0,upper=width,subdivisions=1000L)
   pod<-pod$value/width

   coverage = (sweepWidth*sp)/(width)
   coverages<-c(coverages,coverage)
   pods<-c(pods,pod)
}
dev.off()

print("Normal LRC, parallel sweeps done")

HTMLInsertGraph(graphpcs2,file=output,Caption="Normal lateral range curve exp(-|x^2|), showing multiple parallel sweeps.  One effective sweep width shown in tan.  Cumulative probability in blue.  Each separate sweep in black.")

coveragesSecond <- coverages
podsSecond <- pods

graphpcs3<-"coveragecurves3.png"
png(file.path(outputDirectory,graphpcs3),width=plotW,height=plotH)
par(cex=fontScaling);
x<-c(0,10)
y<-c(0,1)
plot(x,y,type="p",xlim=c(0,10),ylim=c(0,1),main="Normal exp(-|x^10|)*0.9",xlab="Distance",ylab="POD",col="white")

coverages<-c()
pods<-c()
sweepsPerWidth<-c(1,2,3,4,5,6,7,8,9,10)
width<-10
f<-function(x) { exp(-abs(x^10))*0.9 }
area<-integrate(f,-100,100)
sweepWidth<-area$value
f<-function(x,s) { exp(-abs((x-s)^10)) * 0.9 }
dispCov<-0
for (sp in sweepsPerWidth) {
   spacing = width / sp
   offsets<-seq(-1,15)
   offsets<-offsets * spacing
   f.env = new.env()
   f.env$offsets<-offsets
   f.env$f<-f
   if (sp==5) {
      xit<-seq(0,width,by=0.1)
      lines(xit,fToIntegrate(xit),col="blue")
      for (offset in offsets) {
         lines(xit,f(xit,offset))
      }
      polygon(c(spacing-(sweepWidth/2),spacing-(sweepWidth/2),spacing+(sweepWidth/2),spacing+(sweepWidth/2)),c(0,1,1,0),col=adjustcolor("wheat",alpha.f=0.3),border="wheat3")
      dispCov<-(sweepWidth*sp)/width
   }
   pod<-integrate(fToIntegrate,lower=0,upper=width,subdivisions=1000L)
   pod<-pod$value/width

   coverage = (sweepWidth*sp)/(width)
   coverages<-c(coverages,coverage)
   pods<-c(pods,pod)
}
dev.off()

HTMLInsertGraph(graphpcs3,file=output,Caption=paste("Normal lateral range curve exp(-|x^10|)*0.9, showing multiple parallel sweeps at a coverage of ",dispCov,".  One effective sweep width shown in tan.  Cumulative probability in blue.  Each separate sweep in black."))

coveragesThird <- coverages
podsThird <- pods



graphpcs<-"coveragepodsim.png"
png(file.path(outputDirectory,graphpcs),width=plotW,height=plotH)
par(cex=fontScaling);
par(cex=1,cex.lab=fontScaling,cex.axis=fontScaling,mar=c(6, 6, 2, 2)+0.1,lwd=2)
# plot(coveragesFirst,podsFirst,type="l",col="green")
plot(coveragesFirst,podsFirst,type="l",col="green",xlim=c(0,2.0),ylim=c(0,1),xlab="Coverage",ylab="POD")
lines(coverage_exp,pod_exp,col="blue")
lines(coverage_pb,pod_pb,col="black")
lines(coveragesFirst,podsFirst,type="l",col="green",lty="dashed",xlim=c(0,2.0),ylim=c(0,1))
lines(coveragesSecond,podsSecond,type="l",col="red",xlim=c(0,2.0),ylim=c(0,1))
lines(coveragesThird,podsThird,type="l",col="darkred",xlim=c(0,2.0),ylim=c(0,1))
lines(coveragesIter,podsIter,type="l",col="red",lty="dashed",xlim=c(0,2.0),ylim=c(0,1))
lines(coveragesIterPB,podsIterPB,type="l",col="black",lty="dotted",xlim=c(0,2.0),ylim=c(0,1))


coverages<-c()
pods<-c()
sweepsPerWidth<-c(1,2,3,4,5,6,7,8,9,10)
width<-10
f<-function(x) { ifelse( x>=1 | x<=-1 ,0,1) }
area<-integrate(f,-100,100)
sweepWidth<-area$value
f<-function(x,s) { ifelse( x>=1+s | x<=-1+s ,0,1) }
for (sp in sweepsPerWidth) {
   spacing = width / sp
   offsets<-seq(-1,15)
   offsets<-offsets * spacing
   f.env = new.env()
   f.env$offsets<-offsets
   f.env$f<-f
   pod<-integrate(fToIntegrate,lower=0,upper=width,subdivisions=1000L)
   pod<-pod$value/width

   coverage = (sweepWidth*sp)/(width)
   coverages<-c(coverages,coverage)
   pods<-c(pods,pod)
}
lines(coverages,pods,type="l",col="darkgoldenrod4",xlim=c(0,2.0),ylim=c(0,1))

dev.off();

print("Comparison plot done")

rframe <- data.frame(coveragesFirst, podsFirst, coveragesSecond, podsSecond, coveragesThird, podsThird, coveragesIter, podsIter, coveragesIterPB, podsIterPB)
rframe41 <- data.frame(coverage_exp, pod_exp, coverage_pb, pod_pb)
save(rframe,file=file.path(outputDirectory,paste('coveagepodsim_',iterations,'.data',sep='')))
save(rframe41,file=file.path(outputDirectory,paste('coveagepod41sim_',iterations,'.data',sep='')))

HTML("<p>Now let's consider the effect of how sweeps are run against probability of detection.</p>")

HTML("<p>We can have sweeps run parallel in a search segment.  If the parallel sweeps are spaced out one effective sweep width apart (center to center), there is a coverage of 1.  If sweeps are spaced twice the effective sweep width apart (center to center), then there is a coverage of 1/2.  With a perfect broom (definite range) lateral range curve, the relationship between coverage and probability of detection is obvious when sweeps are parallel.  At a coverage of 1, there is a POD of 1.  Placing sweeps closer together to increase the coverage above q1 can't increase the POD above 1.  Placing sweeps further apart linearly decreases the POD.  At a coverage of 1/2, the POD is 50%.  At a coverage of 1/4, the POD is 25%.</p>")

HTML("<p>Things get more interesting if we take that same perfect broom (definite range) lateral range curve and place the sweeps randomly in the search segment.</p>")

HTMLInsertGraph(graphrsf,file=output,Caption="A set of sweeps placed randomly in the search area.")


HTML("<p>Since we are placing the sweeps randomly in the search area, let's repeat multiple times.  Let's take a sensor with a perfect broom (defninite range) lateral range curve placed in increasing numbers of short sweeps in a search area to assess the relationship between coverage and POD.  Each sweep is placed independently of the others (so they can overlap).  Each sweep is short relative to the entire search segment, but each sweep is longer than the definite range of detection.   To get a picture of the relationship between coverage and POD, we'll place 25 sweeps in the search area, then calculate the coverage, excluding from consideration any part of a sweep that falls outside the search area.  We'll then randomly (and independently of each other and the sweeps) place a large number of targets in the search area.  We'll then calculate the POD from the proportion of targets which fall into sweeps and which fall outside any sweep.  This coverage (portion of area swept/portion of area not swept) and POD (proportion of targets in sweeps/proporation of targets not in sweeps) will plot as one point on a POD-coverage graph.   We'll then iterate a few times with 25 sweeps (creating a cluster of a few points).  We'll then repeat, increaseing in sets of 25 sweeps up to 500 sweeps.  This gives a set of points on a coverage - POD graph. </p>")

HTMLInsertGraph(graphrsrsdf,file=output,Caption="Relationship between coverage and POD for simulations of randomly placed definite range sweeps.")

HTML("<p>The points in the graph above for random simulations of detection fall along the line of the exponential detection function - POD (detection) increases as a function of e raised to the negative absolute value of the coverage  p=e^(-|coverage|).    Parallel sweeps with a perfect broom (definite range lateral range function) with a coverage of 1 result in a POD of 1, and increasing coverage can't increase the POD.  Place those sweeps randomly in the search area, and this random placement introduces inefficency in the form of random overlaps - coverage is a measure of the search effort (that is the total area searched), not just the area searce at least once.   </p>")


HTML("<p>Randomly placing relatively short sweeps, each independent of the others, in a search segment, results in coverage gradually increasing - and the more of the area gets covered, the harder it is to randomly place a sweep on an area that hasn't been swept yet.  This results in an exponential relationship between the coverage and the POD.  This exponential detection function is not to be confused with the exponential lateral range curve, even though the two curves have the same formula.  The exponential detection function describes the relationship between coverage and probability of detection for randomly placed sweeps (for any lateral range curve).  The exponential lateral range curve describes visual detection dropping off exponentially with distance from the searcher.</p>")

HTML("<p>So, for a perfect broom (that is, definite range detection function), there are two end members for the detection function (coverage-POD relationship).  If the perfect broom sweeps in neat parallel lines, there is a neat linear relationship between coverage and POD from 0 to 1.  In this range, the coverage equals the POD.  For a coverage over one, POD remains at one.   However, if the perfect broom sweeps are random and independent of each other, the coverage-POD relationship follows the exponential detection function.  Under randomly placed perfect broom sweeps, a coverage of 1 provides a POD of 0.63 (that is, 63%), not one.  With the exponential detection function, low coverages are inefective, and high coverages are inefficeient (as the perfect broom sweeps randomly overlap each other more and more of the time with little gain in added searched area).  The middle ground of coverages around 1 is efficient search.</p>")

HTMLInsertGraph(graphpc,file=output,Caption="Coverage POD relationships (detection functions) the and perfect broom (definite range) lateral range curve under parallel or random sweeps.")



HTML("<p>If we examine the detection functions for a variety of lateral range curves in parallel sweeps, these detection functions end up as curves falling between the detection function for the perfect broom in parallel sweeps and the perfect broom in random sweeps (the exponential detection function).  The exponential detection function can be thought of as a conservative detection function which accounts for random variability and navigation errors in parallel sweeps for any lateral range curve.</p>")

HTMLInsertGraph(graphpcs,file=output,Caption="Coverage POD relationships by simulation for various lateral range curves under parallel sweeps.")

HTMLEndFile()

print("finished")

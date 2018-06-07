# Diagrams of parallel and random sweeps 
#(coverages and pods are for definite range lateral range curve with range of 2).

#X11()
outputDirectory <- file.path("output_figures")
if (!dir.exists(outputDirectory)){ dir.create(outputDirectory) }

# include text for coverage and pod on plot
textOnPlot <- FALSE

plotW<-1200
plotH<-1200

# tan sweeps on blue background 
#sweptColor <- "goldenrod"
#unsweptColor <- "lightskyblue"

# blue sweeps on tan background
sweptColor <- "lightskyblue3"
unsweptColor <- "khaki"
sweepLineColor <- "gray"
alpha <- 0.3

# Functions for plotting sweeps 

plotRandomSweeps <- function(maxiter,doPlot=TRUE,plotCoverage=TRUE) { 

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
     if (doPlot) { lines(c(x1,x2),c(y1,y2),col=sweepLineColor) }
     slope <- - 1/((diff(c(y1,y2))/diff(c(x1,x2))))
     w <- 2
     xes <- c(xm(x1,w,slope),xp(x1,w,slope),xp(x2,w,slope),xm(x2,w,slope))
     yes <- c(ym(y1,w,slope),yp(y1,w,slope),yp(y2,w,slope),ym(y2,w,slope))
     if (doPlot) { polygon(xes,yes,col=adjustcolor(sweptColor,alpha.f=alpha))  }
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
   if (doPlot && plotCoverage) { title(main=paste("coverage=",round(coverage,digits=2)," pod=",round(p,digits=2))) }
   if (doPlot) { plot(spAreaNotSearched,add=TRUE,col=unsweptColor) }

   return( pairlist(p=p,coverage=coverage,podObs=pod) )
}

plotParallelSweeps <- function(totalSweeps,doPlot=TRUE, plotCoverage=TRUE) { 

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
     if (doPlot) { lines(c(x1+searchWidth,x2-searchWidth),c(y1,y2),col=sweepLineColor) }
     slope <- - 1
     w <- 2
     xes <- c(x1,x1,x2,x2)
     yes <- c(y1,y2,y2,y1)
     if (doPlot) { polygon(xes,yes,col=adjustcolor(sweptColor,alpha.f=alpha))  }
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
   if (doPlot && plotCoverage) { title(main=paste("coverage=",round(coverage,digits=2)," pod=",round(p,digits=2))) }
   if (doPlot && !is.null(spAreaNotSearched)) { plot(spAreaNotSearched,add=TRUE,col=unsweptColor) }

   return( pairlist(p=p,coverage=coverage,podObs=pod) )
}

print("setup done")

# plot diagrams of random sweeps

# aproximating coverage of .25
graphrsf<-"randomsweepfrags.png"
png(file.path(outputDirectory,graphrsf),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)
plotRandomSweeps(50,TRUE,textOnPlot)
dev.off();

print("plotted first random sweep")

# aproximating coverage of .5
graphrsf4<-"randomsweepfrags_half.png"
png(file.path(outputDirectory,graphrsf4),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)
plotRandomSweeps(100,TRUE,textOnPlot)
dev.off();

# aproximating coverage of 1
graphrsf2<-"randomsweepfrags_1.png"
png(file.path(outputDirectory,graphrsf2),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)
plotRandomSweeps(190,TRUE,textOnPlot)
dev.off();

# aproximating coverage of 2
graphrsf3<-"randomsweepfrags_2.png"
png(file.path(outputDirectory,graphrsf3),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)
plotRandomSweeps(500,TRUE,textOnPlot)
dev.off();

print("plotted random sweeps")

# plot diagrams of parallel sweeps

graphps<-"parallelsweep.png"
png(file.path(outputDirectory,graphps),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)
plotParallelSweeps(8,TRUE,textOnPlot)
dev.off();


graphps2<-"parallelsweep_half.png"
png(file.path(outputDirectory,graphps2),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)
plotParallelSweeps(13,TRUE,textOnPlot)
dev.off();

graphps3<-"parallelsweep_1.png"
png(file.path(outputDirectory,graphps3),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)
plotParallelSweeps(25,TRUE,textOnPlot)
dev.off();

graphps4<-"parallelsweep_2.png"
png(file.path(outputDirectory,graphps4),width=squarePlotW,height=squarePlotH)
par(cex.main=fontScaling)
plotParallelSweeps(51,TRUE,textOnPlot)
dev.off();


print("done")

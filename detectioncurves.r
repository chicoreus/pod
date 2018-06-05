require("R2HTML")# HTML output

outputDirectory <- file.path("output_detectioncurves")
output <- HTMLInitFile(outputDirectory,filename="pod_output")
HTML.title("Modeling POD");

HTML("This is output from the run of an R program that examines the relationship between detection functions, coverage, and POD.")

#  Distributions
graph1<-"distributioncurves.png"
png(file.path(outputDirectory,graph1),width=1200,height=1200)
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

HTMLInsertGraph(graph1,file=output,Caption="Comparison of 4 detection functions")

# perfect broom
graphbroom<-"broom.png"
png(file.path(outputDirectory,graphbroom),width=1200,height=1200)
f<-function(x) { ifelse( x>=1 | x<=-1 ,0,1) }
y<-f(x)
plot(x,y,type="l",main="Perfect Broom")
area<-integrate(f, -100,100)
lines(c(area$value/2,area$value/2),c(0,1),col="blue",mfg=c(2,2))
lines(c(-area$value/2,-area$value/2),c(0,1),col="blue")
dev.off();

HTMLInsertGraph(graphbroom,file=output,Caption="Perfect broom detection function")

# pod/coverage
graphpc<-"coveragepod.png"
png(file.path(outputDirectory,graphpc),width=1200,height=1200)
#perfect broom
f<-function(x) { ifelse( x>=1 | x<=-1 ,0,1) }
coverage <- seq(0,4,by=0.1)
pod <- c()
for (i in coverage) { 
  area<-integrate(f, 0,i)
  pod<-c(pod,area$value)
}
plot(coverage,pod,type="l")
coverage_pb<-coverage
pod_pb<-pod
# exponential 
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

HTMLInsertGraph(graphpc,file=output,Caption="Coverage POD relationships for exponetnials and perfect broom detection fuctions")

graphpcs1<-"coveragecurves.png"
png(file.path(outputDirectory,graphpcs1),width=1200,height=1200)

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
plot(x,y,type="p",xlim=c(0,10),ylim=c(0,1),main="Exponential detection function")

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
   hits=0
   tries=10000
   for (i in 1:tries) { 
      x <- runif(1,min=0,max=width)
      thresholds <- lp(f,x,spacing)

if (sp==1) { 
xes<-c()
for(threshold in thresholds) { 
   xes<-c(xes,x)
}
points(xes,thresholds)
}
if (sp==5) { 
xes<-c()
for(threshold in thresholds) { 
   xes<-c(xes,x)
}
points(xes,thresholds,col="blue")
}

#      detect<-sum(cp.quadratic(thresholds)[-1]) 
#     same result as 
      detect = pUnionAnyOneOf(thresholds)
      p<-runif(1)
      if (p<detect) { hits<-hits+1 } 
   } 
   coverage = (sweepWidth*sp)/(width)
   pod<-(hits/tries)
   coverages<-c(coverages,coverage)
   pods<-c(pods,pod)
} 
dev.off()

coveragesFirst <- coverages
podsFirst <- pods


#coverages<-c()
#pods<-c()
#widths<-c(1.5,2,2.5,3,3.5,4,5,6,7,8,9,10,11,12,14,16,20,24,28,36,46,66)
#f<-function(x) {exp(-abs(x^2))}
#area<-integrate(f,-100,100)
#halfSweepWidth<-area$value/2
#for (width in widths) {
#   hits=0
#   misses=0
#   spacing = width/2
#   f1<-function(x) {exp(-abs((x-spacing)^2))}
#   f2<-function(x) {exp(-abs(((x-(spacing*2)^2))))}
#   f3<-function(x) {exp(-abs(((x-(spacing*3)^2))))}
#   f4<-function(x) {exp(-abs(((x-(spacing*-1)^2))))}
#   for (i in 1:10000) {
#      x <- runif(1,0,width)
#      detect <- f(x)
#      detect1 <- f1(x)
#      detect2 <- f2(x)
#      detect3 <- f3(x)
#      detect4 <- f4(x)
#      p<-runif(1,0,1)
#      if (p<detect) { hits<-hits+1 } else {
#         if (p<detect1) { hits<-hits+1 } else {
#           if (p<detect2) { hits<-hits+1 } else {
#           if (p<detect3) { hits<-hits+1 } else {
#           if (p<detect4) { hits<-hits+1 } else {
#             misses<-misses+1
#           }
#           }
#           }
#         }
#      }
#   }
#   coverage = (halfSweepWidth*4)/(1*width)
#   pod<- hits/(hits+misses)
#   coverages<-c(coverages,coverage)
#   pods<-c(pods,pod)
#}
#lines(coverages,pods,type="l",col="red")

graphpcs2<-"coveragecurves2.png"
png(file.path(outputDirectory,graphpcs2),width=1200,height=1200)
x<-c(0,10)
y<-c(0,1)
plot(x,y,type="p",xlim=c(0,10),ylim=c(0,1),main="Normal exp(-|x^2|)")

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
   hits=0
   tries=10000
   for (i in 1:tries) {
      x <- runif(1,min=0,max=width)
      thresholds <- lp(f,x,spacing)

if (sp==1) {
xes<-c()
for(threshold in thresholds) {
   xes<-c(xes,x)
}
points(xes,thresholds)
}
if (sp==5) {
xes<-c()
for(threshold in thresholds) {
   xes<-c(xes,x)
}
points(xes,thresholds,col="blue")
}

#      detect<-sum(cp.quadratic(thresholds)[-1]) 
      detect = pUnionAnyOneOf(thresholds)
      p<-runif(1)
      if (p<detect) { hits<-hits+1 }
   }
   coverage = (sweepWidth*sp)/(width)
   pod<-(hits/tries)
   coverages<-c(coverages,coverage)
   pods<-c(pods,pod)
}
dev.off()

coveragesSecond <- coverages
podsSecond <- pods


graphpcs3<-"coveragecurves3.png"
png(file.path(outputDirectory,graphpcs3),width=1200,height=1200)
x<-c(0,10)
y<-c(0,1)
plot(x,y,type="p",xlim=c(0,10),ylim=c(0,1),main="Normal exp(-|x^10|)*0.9")

coverages<-c()
pods<-c()
sweepsPerWidth<-c(1,2,3,4,5,6,7,8,9,10)
width<-10
f<-function(x) { exp(-abs(x^10))*0.9 }
area<-integrate(f,-100,100)
sweepWidth<-area$value
f<-function(x,s) { exp(-abs((x-s)^10)) * 0.9 }
for (sp in sweepsPerWidth) {
   spacing = width / sp
   hits=0
   tries=10000
   for (i in 1:tries) {
      x <- runif(1,min=0,max=width)
      thresholds <- lp(f,x,spacing)

if (sp==1) {
xes<-c()
for(threshold in thresholds) {
   xes<-c(xes,x)
}
points(xes,thresholds)
}
if (sp==5) {
xes<-c()
for(threshold in thresholds) {
   xes<-c(xes,x)
}
points(xes,thresholds,col="blue")
}

#      detect<-sum(cp.quadratic(thresholds)[-1]) 
      detect = pUnionAnyOneOf(thresholds)
      p<-runif(1)
      if (p<detect) { hits<-hits+1 }
   }
   coverage = (sweepWidth*sp)/(width)
   pod<-(hits/tries)
   coverages<-c(coverages,coverage)
   pods<-c(pods,pod)
}
dev.off()

coveragesThird <- coverages
podsThird <- pods



graphpcs<-"coveragepodsim.png"
png(file.path(outputDirectory,graphpcs),width=1200,height=1200)
plot(coveragesFirst,podsFirst,type="l",col="green",xlim=c(0,2.0),ylim=c(0,1))
lines(coverage_exp,pod_exp,col="blue")
lines(coverage_pb,pod_pb,col="black")
lines(coveragesSecond,podsSecond,type="l",col="red",xlim=c(0,2.0),ylim=c(0,1))
lines(coveragesThird,podsThird,type="l",col="red",xlim=c(0,2.0),ylim=c(0,1))

#
#coverages<-c()
#pods<-c()
#widths<-c(1.5,2,2.5,3,3.5,4,5,6,7,8,9,10,11,12,14,16,20,24,28,36,46,66)
#f<-function(x) { ifelse( x>=1 | x<=-1 ,0,1) }
#area<-integrate(f,-100,100)
#halfSweepWidth<-area$value/2
#for (width in widths) {
#   hits=0
#   misses=0
#   spacing = width/2
#   f1<-function(x) { ifelse( x>=1+spacing | x<=-1+spacing ,0,1) }
#   f2<-function(x) { ifelse( x>=1+spacing*2 | x<=-1+2*spacing ,0,1) }
#   f3<-function(x) { ifelse( x>=1+spacing*3 | x<=-1+3*spacing ,0,1) }
#   f4<-function(x) { ifelse( x>=1-spacing | x<=-1-spacing ,0,1) }
#   for (i in 1:10000) {
#      x <- runif(1,0,width)
#      detect <- f(x)
#      detect1 <- f1(x)
#      detect2 <- f2(x)
#      detect3 <- f3(x)
#      detect4 <- f4(x)
#      p<-runif(1,0,1)
#      if (p<detect) { hits<-hits+1 } else {
#         if (p<detect1) { hits<-hits+1 } else {
#           if (p<detect2) { hits<-hits+1 } else {
#           if (p<detect3) { hits<-hits+1 } else {
#           if (p<detect4) { hits<-hits+1 } else {
#             misses<-misses+1
#           }
#           }
#           }
#         }
#      }
#   }
#   coverage = (halfSweepWidth*4)/(1*width)
#   pod<- hits/(hits+misses)
#   coverages<-c(coverages,coverage)
#   pods<-c(pods,pod)
#}
#lines(coverages,pods,type="l",col="blue")

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
   hits=0
   tries=10000
   for (i in 1:tries) { 
      x <- runif(1,min=0,max=width)
      thresholds <- lp(f,x,spacing)
#      detect<-sum(cp.quadratic(thresholds)[-1]) 
      detect = pUnionAnyOneOf(thresholds)
      p<-runif(1)
      if (p<detect) { hits<-hits+1 } 
   } 
   coverage = (sweepWidth*sp)/(width)
   pod<-(hits/tries)
   coverages<-c(coverages,coverage)
   pods<-c(pods,pod)
} 
lines(coverages,pods,type="l",col="green",xlim=c(0,2.0),ylim=c(0,1))

dev.off();
HTMLInsertGraph(graphpcs,file=output,Caption="Coverage POD relationships by simulation for various detection fuctions")

HTMLEndFile()

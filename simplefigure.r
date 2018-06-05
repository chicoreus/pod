# A couple of simple figures.

#X11()
outputDirectory <- file.path("output")
if (!dir.exists(outputDirectory)){ dir.create(outputDirectory) }

plotW<-1200
plotH<-1200

# Comparisons of existing data/extrapolated curves for canine POD.
# From Graham's original canine POD work, from the NASAR MLPI text, 
# and from Chiacchia et al., 2015.

graph2<-"podcomparison_simple.png"
png(file.path(outputDirectory,graph2),width=plotW,height=plotH)

# Setup graphical parameters
par(mfrow=c(1,1))
par(cex=3)
par(lwd=2)
par(bty='l')

# NASAR 
n100m<-c(10,26,39)
n50m<-c(18,46,63)
stability<-c(1,3,4)
nasar<-data.frame(stability,n100m,n50m)
plot(c(1,3,4),nasar[,2],xlim=c(1,7.7),ylim=c(0,100),xaxt="n",type='l',xlab='Atmospheric Stability Category',ylab='POD (%)')
axis(1,1:6,LETTERS[1:6])
lines(c(1,3,4),nasar[,3],xlim=c(1,6),ylim=c(0,100),xaxt="n",lty=5)
text(c(4,4),c(39,63),labels=c('MLPI 100m','MLPI 50m'),pos=4)
text(c(5,5),c(10,20),labels=c('A=Least Stable','F=Most Stable'),pos=4)

# Graham
stability <- c(1,2,3,4,5,6)
g100m <- c(5,10,35,80,90,95)
g50m <- c(50,55,67,90,95,97)
graham<-data.frame(stability,g100m,g50m)
lines(c(1,2,3,4,5,6),graham[,2],xlim=c(1,6),ylim=c(0,100),xaxt="n",col='blue')
lines(c(1,2,3,4,5,6),graham[,3],xlim=c(1,6),ylim=c(0,100),xaxt="n",col='blue',lty=5)
points(c(1,6),c(graham[1,2],graham[6,2]),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='blue')
text(c(6),c(95),labels=c('Graham 100m'),pos=4,col='blue')

# Chiacchia
lines(c(2,3,4),c(55,63,70),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red')
lines(c(2,3,4),c(80,86,91),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red',lty=5)
points(c(2,3,4),c(80,86,91),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red')
text(c(4,4),c(70,91),labels=c('Chiacchia 100m','Chiacchia 50m'),pos=4,col='red')

dev.off()

# Graph of an exponential lateral range curve and an exponential lateral
# range curve with an exponent of 10 ( e^(-|x|^10) ) but multiplied by 0.9
graph1<-"exp_exp10_LRC_curves.png"
png(file.path(outputDirectory,graph1),width=plotW,height=plotH)
par(mfrow=c(2,1))
par(lwd=2)
par(cex=1.2)
par(mgp=c(5,1,0))
# Exponential lateral range curve 
x<-seq(-3,3,0.01)
y<-exp(-abs(x))
plot(x,y,type="l", ylab=list("POD", cex=3), xlab=list("", cex=.1))
aoc<-sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
lines(c(aoc/2,aoc/2),c(0,1),col='blue')
lines(c(-aoc/2,-aoc/2),c(0,1),col='blue')
# Approximating perfect broom, but peak at 0.9 
y<-(exp(-abs(x)^10)) * 0.9
plot(x,y,type="l",ylim=c(0,1), ylab=list("POD",cex=3), xlab=NULL)
aoc<-sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
lines(c(aoc/2,aoc/2),c(0,1),col='blue')
lines(c(-aoc/2,-aoc/2),c(0,1),col='blue')
dev.off();


# Same graph, but with areas shaded.
graph2<-"exp_exp10_LRC_curvesshaded.png"
png(file.path(outputDirectory,graph2),width=plotW,height=plotH)
par(mfrow=c(2,1))
par(lwd=2)
par(cex=1.2)
par(cex.axis=2.2)
# Exponential lateral range curve 
x<-seq(-3,3,0.01)
y<-exp(-abs(x))
plot(x,y,type="h", ylab="POD", xlab=NULL,col="grey95",ann=FALSE,ylim=c(0,1))
lines(x,y,type="l",col="black")
aoc<-sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
lines(c(aoc/2,aoc/2),c(0,1),col='blue')
lines(c(-aoc/2,-aoc/2),c(0,1),col='blue')
arx<-c(aoc/2,aoc/2,-aoc/2,-aoc/2,aoc/2)
ary<-c(0,1,1,0,0)
polygon(arx,ary,col=adjustcolor("wheat",alpha.f=0.3),border=NA)
arrows(aoc/2, 0.6, -aoc/2 ,0.6, col="blue",code=3)
# Approximating perfect broom, but peak at 0.9 
y<-(exp(-abs(x)^10)) * 0.9
plot(x,y,type="h",ylim=c(0,1), ylab="POD", xlab=NULL,col="grey95",ann=FALSE)
lines(x,y,type="l",col="black")
aoc<-sum(diff(x) * (head(y,-1)+tail(y,-1)))/2
lines(c(aoc/2,aoc/2),c(0,1),col='blue')
lines(c(-aoc/2,-aoc/2),c(0,1),col='blue')
arx<-c(aoc/2,aoc/2,-aoc/2,-aoc/2,aoc/2)
ary<-c(0,1,1,0,0)
polygon(arx,ary,col=adjustcolor("wheat",alpha.f=0.3),border=NA)
arrows(aoc/2, 0.6, -aoc/2 ,0.6, col="blue",code=3)
dev.off();


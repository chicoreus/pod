# A couple of simple figures.

#X11()
outputDirectory <- file.path("output")
if (!dir.exists(outputDirectory)){ dir.create(outputDirectory) }

plotW<-1200
plotH<-1200

# obtain estimates for POD for sweep width from Chiacchia's 2015 study.
# Figure 4 shows a regression between air stability classes 1-4 (equivalent to A-D)
# and effective sweep width: ESW estimated at 29*AirStabilityCategory + 28 
# Regression with slope of 29 m/airstabilitycategory and ESW intercept of 28
# Lots of caveats - r^2 is .47, lots of spread in the data.  Limited to a 
# similarly trained small number of dogs in one region.  
# Note: Limited data for category A, excluded from plots below.

airStability <- c(1,2,3,4)
regression <- function(x) { (29*x)+28 } 
ESW <- regression(airStability)
# calculate coverages for these ESWs at given sweep widths
coverage25 <- ESW/25
coverage50 <- ESW/50
coverage100 <- ESW/100
# assuming an exponential detection function, calculate the POD for the coverage: 
expDF <- function(x) { 1-exp(-abs(x))}
pod25 <- expDF(coverage25)
pod50 <- expDF(coverage50)
pod100 <- expDF(coverage100)
# describe in a function that takes an airstability class and a sweepwidth and returns a pod in percent
chia <- function(air,sweepw) { expDF(regression(air)/sweepw) * 100 }


# Comparisons of existing data/extrapolated curves for canine POD.
# From Graham's original canine POD work, from the NASAR MLPI text, 
# and from Chiacchia et al., 2015.

graph<-"podcomparison_simple.png"
png(file.path(outputDirectory,graph),width=plotW,height=plotH)

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



# Comparisons of existing data/extrapolated curves for canine POD.
# From Graham's original canine POD work, from the NASAR MLPI text, 
# and from Chiacchia et al., 2015.
# Including a proposal for values and spacings to use combiniing these.

plotW<-1500

graph3<-"podcomparison_proposal.png"
png(file.path(outputDirectory,graph3),width=plotW,height=plotH)

# Setup graphical parameters
par(mfrow=c(1,1))
par(cex=3)
par(lwd=2)
par(bty='l')

# NASAR ,
n100m<-c(10,26,39)
n50m<-c(18,46,63)
n25m<-c(33,70,86)
stability<-c(1,3,4)
nasar<-data.frame(stability,n100m,n50m,n25m)
plot(c(1,3,4),nasar[,2],xlim=c(1,7.7),ylim=c(0,100),xaxt="n",type='l',xlab='Atmospheric Stability Category',ylab='POD (%)')
axis(1,1:6,LETTERS[1:6])
lines(c(1,3,4),nasar[,3],xlim=c(1,6),ylim=c(0,100),xaxt="n",lty=5)
lines(c(1,3,4),nasar[,4],xlim=c(1,6),ylim=c(0,100),xaxt="n",lty=3)
#text(c(4,4),c(39,63),labels=c('MLPI 100m','MLPI 50m'),pos=4)
text(c(3.25,3.25),c(10,20),labels=c('A=Least Stable','F=Most Stable'),pos=4)

# Graham
stability <- c(1,2,3,4,5,6)
g100m <- c(5,10,35,80,90,95)
g50m <- c(50,55,67,90,95,97)
g25m <- c(75,77,86,95,97,99)
graham<-data.frame(stability,g100m,g50m,g25m)
lines(c(1,2,3,4,5,6),graham[,2],xlim=c(1,6),ylim=c(0,100),xaxt="n",col='blue')
lines(c(1,2,3,4,5,6),graham[,3],xlim=c(1,6),ylim=c(0,100),xaxt="n",col='blue',lty=5)
lines(c(1,2,3,4,5,6),graham[,4],xlim=c(1,6),ylim=c(0,100),xaxt="n",col='blue',lty=3)
points(c(1,6),c(graham[1,2],graham[6,2]),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='blue')
#text(c(6),c(95),labels=c('Graham 100m'),pos=4,col='blue')


# Chiacchia
lines(c(2,3,4),c(chia(2,100),chia(3,100),chia(4,100)),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red')
lines(c(2,3,4),c(chia(2,50),chia(3,50),chia(4,50)),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red',lty=5)
lines(c(2,3,4),c(chia(2,25),chia(3,25),chia(4,25)),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red',lty=3)
# We'll treet just one of the Chiacchia numbers as data, though that isn't quite true there is a lot of
# estimatioo between the sweep width data and these numbers, and they have caveats on their study,
# so use the overall crossover point of 50 meters for the full data set as the place to put data markers
points(c(2,3,4),c(chia(2,50),chia(3,50),chia(4,50)),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red')
#> pod25
#[1] 0.8977158 0.9679353 0.9899482 0.9968489
#> pod50
#[1] 0.6801810 0.8209339 0.8997412 0.9438652
#> pod100
#[1] 0.4344746 0.5768379 0.6833632 0.7630722
# estimations in these numbers seem a bit off, not clear why.
#lines(c(2,3,4),c(55,63,70),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red')
#lines(c(2,3,4),c(80,86,91),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red',lty=5)
#lines(c(2,3,4),c(96,98,99),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red',lty=3)
#points(c(2,3,4),c(80,86,91),xlim=c(1,6),ylim=c(0,100),xaxt="n",col='red')
#text(c(4,4),c(70,91),labels=c('Chiacchia 100m','Chiacchia 50m'),pos=4,col='red')

points(c(1),c(55),xlim=c(1,6),ylim=c(0,100),xaxt="n",bg='yellow3',pch=22)
text(c(1),c(55),labels=c('25m'),pos=4,col='yellow3')
points(c(2),c(70),xlim=c(1,6),ylim=c(0,100),xaxt="n",bg='yellow3',pch=22)
text(c(2),c(70),labels=c('75m'),pos=4,col='yellow3')
points(c(3),c(63),xlim=c(1,6),ylim=c(0,100),xaxt="n",bg='brown',pch=22)
text(c(3),c(63),labels=c('100m'),pos=4,col='brown')
points(c(4),c(70),xlim=c(1,6),ylim=c(0,100),xaxt="n",bg='green',pch=22)
text(c(4),c(70),labels=c('100m'),pos=4,col='green')
points(c(5),c(86),xlim=c(1,6),ylim=c(0,100),xaxt="n",bg='green',pch=22)
text(c(5),c(86),labels=c('100m'),pos=4,col='green')
points(c(6),c(86),xlim=c(1,6),ylim=c(0,100),xaxt="n",bg='green',pch=22)
text(c(6),c(86),labels=c('100m'),pos=4,col='green')

legend("bottomright", pch = c(NA,NA,NA,  NA,NA,NA,1,  NA,NA,NA,1, 22), 
                      lty=c(1,5,3,  1,5,3,NA,  1,3,5,NA, NA), 
                      col = c("black", "black", "black",  "blue", "blue", "blue", "blue", "red","red","red","red", "black"), 
        legend = c("NASAR MLPI 100m", "NASAR MLPI 50m", "NASAR MLPI 25m", 
                   "Graham 100m", "Graham 50m", "Graham 25m", "Graham Data", 
                   "Chiacchia 100m","Chiacchia 50m","Chiacchia 25m","Chiacchia Data",   "Proposed Targets"))

dev.off()



pkgname <- "normalp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('normalp')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dnormp")
### * dnormp

flush(stderr()); flush(stdout())

### Name: dnormp
### Title: Density function of an exponential power distribution
### Aliases: dnormp
### Keywords: distribution

### ** Examples

## Compute the density for a vector x with mu=0, sigmap=1 and p=1.5
## At the end we have the graph of the exponential power distribution 
## density function with p=1.5
x <- c(-1, 1)
f <- dnormp(x, p=1.5)
print(f)
plot(function(x) dnormp(x, p=1.5) , -4, 4,
          main = "Exponential power distribution density function (p=1.5)", ylab="f(x)")



cleanEx()
nameEx("estimatep")
### * estimatep

flush(stderr()); flush(stdout())

### Name: estimatep
### Title: Estimation of p
### Aliases: estimatep
### Keywords: univar

### ** Examples

x<-rnormp(300,mu=1,sigmap=2,p=4)
p<-estimatep(x,mu=1,p=2)
p



cleanEx()
nameEx("graphnp")
### * graphnp

flush(stderr()); flush(stdout())

### Name: graphnp
### Title: Plot of exponential power distributions
### Aliases: graphnp
### Keywords: aplot

### ** Examples

## Plot four different curves with p=1,2,3,4
## and 50 (it will plot an uniform distribution)
graphnp(c(1:4,50))



cleanEx()
nameEx("kurtosis")
### * kurtosis

flush(stderr()); flush(stdout())

### Name: kurtosis
### Title: Indices of kurtosis
### Aliases: kurtosis
### Keywords: univar

### ** Examples

kurtosis(p=2)
x<-rnormp(50,mu=0,sigmap=2,p=1.5)
kurtosis(x,p=2)



cleanEx()
nameEx("lmp")
### * lmp

flush(stderr()); flush(stdout())

### Name: lmp
### Title: Fitted linear model with exponential power distribution errors
### Aliases: lmp lmp.default
### Keywords: regression

### ** Examples

e<-rnormp(n=100,mu=0,sigmap=4,p=3,method="d")
x<-runif(100)
y<-0.5+2*x+e
lmp(y~x)



cleanEx()
nameEx("paramp")
### * paramp

flush(stderr()); flush(stdout())

### Name: paramp
### Title: Estimation of location and scale parameters
### Aliases: paramp paramp.default print.paramp
### Keywords: univar

### ** Examples

x<-rnormp(1000,2,3,4.2)
paramp(x)



cleanEx()
nameEx("plot.lmp")
### * plot.lmp

flush(stderr()); flush(stdout())

### Name: plot.lmp
### Title: Diagnostic plots for a lmp object
### Aliases: plot.lmp
### Keywords: hplot

### ** Examples

x<-1:20
z<-runif(20)
e<-rnormp(20,mu=0,sigmap=1,p=3)
y<-0.5+x+z+e
lmp.res<-lmp(y~x+z)
plot(lmp.res)



cleanEx()
nameEx("plot.simul.lmp")
### * plot.simul.lmp

flush(stderr()); flush(stdout())

### Name: plot.simul.lmp
### Title: Plots of the results of a simulation plan on a linear regression
###   model
### Aliases: plot.simul.lmp
### Keywords: hplot

### ** Examples

sim<-simul.lmp(n=10,m=50,q=1,data=1.5,int=0,sigmap=1,p=3.5)
plot(sim)



cleanEx()
nameEx("plot.simul.mp")
### * plot.simul.mp

flush(stderr()); flush(stdout())

### Name: plot.simul.mp
### Title: Plots of the results of a simulation plan on the parameters of
###   an exponential power distribution
### Aliases: plot.simul.mp
### Keywords: hplot

### ** Examples

## The histograms of all the computed estimates
a<-simul.mp(100,50,mu=0,sigmap=1,p=3)
plot(a)



cleanEx()
nameEx("pnormp")
### * pnormp

flush(stderr()); flush(stdout())

### Name: pnormp
### Title: Probability function of an exponential power distribution
### Aliases: pnormp
### Keywords: distribution

### ** Examples

## Compute the distribution function for a vector x with mu=0, sigmap=1 and p=1.5
## At the end we have the graph of the exponential power distribution function with p=1.5.
x <- c(-1, 1)
pr <- pnormp(x, p=1.5)
print(pr)
plot(function(x) pnormp(x, p=1.5), -4, 4,
          main = "Exponential Power Distribution Function (p=1.5)", ylab="F(x)")



cleanEx()
nameEx("qnormp")
### * qnormp

flush(stderr()); flush(stdout())

### Name: qnormp
### Title: Quantiles of an exponential power distribution
### Aliases: qnormp
### Keywords: distribution

### ** Examples

## Compute the quantiles for a vector of probabilities x
## with mu=1, sigmap=2 and p=1.5
x <- 0.3
q <- qnormp(x, 1, 2, 1.5)
q



cleanEx()
nameEx("qqnormp")
### * qqnormp

flush(stderr()); flush(stdout())

### Name: qqnormp
### Title: Quantile-Quantile plot for an exponential power distribution
### Aliases: qqnormp qqnormp.default qqlinep
### Keywords: hplot

### ** Examples

## Exponential power distribution Q-Q plot for a sample of 100 observations.
e<-rnormp(100,mu=0,sigmap=1,p=3)
qqnormp(e,p=3)
qqlinep(e,p=3)



cleanEx()
nameEx("rnormp")
### * rnormp

flush(stderr()); flush(stdout())

### Name: rnormp
### Title: Pseudo-random numbers from an exponential power distribution
### Aliases: rnormp
### Keywords: distribution

### ** Examples

## Generate a random sample x from an exponential power distribution
## At the end we have the histogram of x
x <- rnormp(1000, 1, 2, 1.5)
hist(x, main="Histogram of the random sample")



cleanEx()
nameEx("simul.lmp")
### * simul.lmp

flush(stderr()); flush(stdout())

### Name: simul.lmp
### Title: Simulation planning for a linear regression model with errors
###   distributed as an exponential power distribution
### Aliases: simul.lmp
### Keywords: regression

### ** Examples

## Simulation of 50 samples of size 10 for a linear regression model with 1 regressor.
simul.lmp(10,50,1,data=1.5,int=1,sigmap=1,p=3,lp=FALSE)



cleanEx()
nameEx("simul.mp")
### * simul.mp

flush(stderr()); flush(stdout())

### Name: simul.mp
### Title: Simulation planning for the parameters of an exponential power
###   distribution
### Aliases: simul.mp print.simul.mp simul.mp.default
### Keywords: univar

### ** Examples

## Simulation plan for 100 samples of size 20, with mu=0, sigmap=1, p=3.
simul.mp(20,100,mu=0,sigmap=1,p=3)



cleanEx()
nameEx("summary.lmp")
### * summary.lmp

flush(stderr()); flush(stdout())

### Name: summary.lmp
### Title: Summarize linear model fits with exponential power distribution
###   errors
### Aliases: summary.lmp print.summary.lmp
### Keywords: regression

### ** Examples

x<-runif(30)
e<-rnormp(30,0,3,1.25)
y<-0.5+x+e
L<-lmp(y~x)
summary(L)



cleanEx()
nameEx("summary.simul.lmp")
### * summary.simul.lmp

flush(stderr()); flush(stdout())

### Name: summary.simul.lmp
### Title: Summarize simulation results on linear regression model
### Aliases: summary.simul.lmp print.summary.simul.lmp
### Keywords: regression

### ** Examples

ris<-simul.lmp(100,20,2,data=c(3,2),int=0,sigmap=1,p=3)
summary(ris)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

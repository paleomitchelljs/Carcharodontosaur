library(geiger)
library(phytools)
library(convevol)

source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\cal3_Functions.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\placeFossil.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\addScale.R")

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")

Input <- readRDS("staytonCC_output.rds")

allTax <- sapply(Input[[1]], function(x) x[,2])
Metrics <- sapply(Input[[1]], function(x) x[,1])
merTax <- sapply(Input[[2]], function(x) x[,2])
abelTax <- sapply(Input[[3]], function(x) x[,2])

pallTax <- sapply(Input[[5]], function(x) x[,2])
pmerTax <- sapply(Input[[6]], function(x) x[,2])
pabelTax <- sapply(Input[[7]], function(x) x[,2])

Plot <- FALSE
if (Plot)	{
pdf("allConvPvals.pdf", height = 5, width = 10)
Breaks <- seq(from=0, to=1, by=0.025)
par(mfrow=c(1,2), bty="l", las=1, mgp=c(1.75, 0.25, 0), tck=-0.01, mar=c(3,4,1,1))
Alpha <- 0.5
Cols <- c(rgb(27/255,158/255,119/255, Alpha), rgb(217/255,95/255,2/255, Alpha), rgb(117/255,112/255,179/255,Alpha))
hist(allTax, xlim=c(0, 1), ylim=c(0,300), main="forelimb & body size", col=Cols[1], border=NA, breaks=Breaks, xlab="p-value", ylab="")
hist(merTax, xlim=c(0, 1), ylim=c(0,300), main="", col=Cols[2], border=NA, add=T, breaks=Breaks, xlab="", ylab="")
hist(abelTax, xlim=c(0, 1), ylim=c(0,300), main="", col=Cols[3], border=NA, add=T, breaks=Breaks, xlab="", ylab="")
abline(v=0.05, lty=3, lwd=2)
Alpha <- 1
Cols <- c(rgb(27/255,158/255,119/255, Alpha), rgb(217/255,95/255,2/255, Alpha), rgb(117/255,112/255,179/255,Alpha))
legend("topright", legend=c("tyranno, abel, carcharo", "Meraxes - Tarbosaurus", "Meraxes - Carnotaurus"), col=Cols, pch=15, bty="n")

Alpha <- 0.5
Cols <- c(rgb(27/255,158/255,119/255, Alpha), rgb(217/255,95/255,2/255, Alpha), rgb(117/255,112/255,179/255,Alpha))
hist(pallTax, xlim=c(0, 1), ylim=c(0,150), col=Cols[1], border=NA, main="pPCA convergence", breaks=Breaks, xlab="p-value", ylab="")
hist(pmerTax, xlim=c(0, 1), ylim=c(0,150), main="", col=Cols[2], border=NA, add=T, breaks=Breaks, xlab="", ylab="")
hist(pabelTax, xlim=c(0, 1), ylim=c(0,150), main="", col=Cols[3], border=NA, add=T, breaks=Breaks, xlab="", ylab="")
abline(v=0.05, lty=3, lwd=2)
Alpha <- 1
Cols <- c(rgb(27/255,158/255,119/255, Alpha), rgb(217/255,95/255,2/255, Alpha), rgb(117/255,112/255,179/255,Alpha))
legend("topright", legend=c("tyranno, abel, carcharo", "Meraxes - Tarbosaurus", "Meraxes - Carnotaurus"), col=Cols, pch=15, bty="n")
dev.off()
}

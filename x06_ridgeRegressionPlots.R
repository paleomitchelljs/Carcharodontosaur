setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")
library(RRphylo)
library(phytools)

source("x00_readData.R")

ACES <- read.csv("rrphylo_aces.csv", row.names=1)
rMat <- read.csv("rrphylo_rates.csv", row.names=1)
Cors <- c(read.csv("rrphylo_cors.csv")[,2])
trees <- read.tree("rrphylo_trees.tre")
analysisTree <- trees[[which(Cors == max(Cors))]]

Arms <- read.csv("carch_measures.csv", stringsAsFactors = F, row.names = 1)
Data <- cbind(apply(Arms[analysisTree$tip.label,2:4], 1, sum), Arms[analysisTree$tip.label,5])

ArmL <- c(Data[,1])
FemurL <- c(Data[,2])
Ratio <- log10(ArmL) - log10(FemurL)

Dat <- c(apply(ACES, 1, weighted.mean, w=Cors), Ratio)
names(Dat) <- c(63:123, names(Ratio))

#rateMap <- contMap(analysisTree, apply(rMat[analysisTree$tip.label,], 1, weighted.mean, w=Cors^2), method="user", anc.states=apply(rMat[rownames(ACES),], 1, weighted.mean, w=Cors^2))
#rateMap <- setMap(rateMap, c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b'))
Map <- contMap(analysisTree, Ratio, method="user", anc.states=apply(ACES, 1, weighted.mean, w=Cors^2), plot=F)

#Map <- setMap(Map, c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'))

# mapping for ratio contmap plot
Map <- setMap(Map, c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b'))
#Map <- setMap(Map, c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

#pdf("armRatio+rateContMap.pdf", height=4.5, width=9)
#par(mfrow=c(1,2), las=1, mgp=c(2, 0.5, 0), tck=-0.01)
#plot(Map, fsize=0.4, sig=2, leg.txt="log Arm - log Femur")
#plot(rateMap, fsize=0.4, sig=2, leg.txt="rate")
#dev.off()

#alpha <- 0.5
#COLS <- c(rgb(166/255,206/255,227/255,alpha), rgb(178/255,223/255,138/255, alpha))
#hist(apply(rMat,1,weighted.mean, w=Cors^2), main="", xlab="", xlim=c(-0.01, 0.01), col=COLS[1])
#Meraxes <- hist(as.numeric(rMat["Meraxes",]), add=T, col=COLS[2])

#plot(log10(Data[analysisTree$tip.label,2]), apply(rMat[analysisTree$tip.label,], 1, weighted.mean, w=Cors^2), pch=16)

#RRtrend <- search.trend(RR=Ridge, y=Ratio, foldername="C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")
#RRconv <- search.conv(RR=Ridge, y=Data)

#RRtrend <- search.trend(RR=Ridge, y=ArmL, x1=FemurX, foldername="C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")

#RidgeFemur <- RRphylo(analysisTree, FemurL)
#FemurX <- c(RidgeFemur$aces[,1], FemurL)
#Ridge <- RRphylo(analysisTree, ArmL, x1=FemurX)
#contMap(analysisTree, ArmL, anc.states=Ridge$aces[rownames(ACES),1])
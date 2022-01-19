setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")
library(RRphylo)
library(phytools)

source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\cal3_Functions.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\placeFossil.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\addScale.R")

Nruns <- 100
Cors <- c()
SSE <- c()
ACES <- matrix(NA, nrow=61, ncol=Nruns)
rMat <- matrix(NA, nrow=123, ncol=Nruns)
trees <- list()
for (i in 1:Nruns)	{
	source("x00_readData.R")
	trees[[i]] <- analysisTree
	Data <- cbind(apply(Arms[analysisTree$tip.label,2:4], 1, sum), Arms[analysisTree$tip.label,5])
	Data2 <- Arms[analysisTree$tip.label,2:5]

	ArmL <- c(Data[,1])
	FemurL <- c(Data[,2])

	Ratio <- log10(ArmL) - log10(FemurL)
	Ridge <- RRphylo(analysisTree, Ratio)

	rMat[,i] <- Ridge$rates
	ACES[,i] <- Ridge$aces
	Cors[i] <- cor(Ridge$predicted.phenotype, Ratio)
	SSE[i] <- sum((Ratio - Ridge$predicted.phenotype)^2)
	cat(i, " ", sum(analysisTree$edge.length), "\n")
}

rownames(ACES) <- rownames(Ridge$aces)
rownames(rMat) <- rownames(Ridge$rates)

write.csv(ACES, "rrphylo_aces.csv")
write.csv(rMat, "rrphylo_rates.csv")
write.csv(Cors, "rrphylo_cors.csv")
class(trees) <- "multiPhylo"
write.tree(trees, "rrphylo_trees.tre")
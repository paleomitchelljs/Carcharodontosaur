library(phytools)

source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\cal3_Functions.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\placeFossil.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\addScale.R")

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")
source("x00_readData.R")

#### pPCA
Data <- log10(Arms[analysisTree$tip.label,2:5])
pPCA <- phyl.pca(analysisTree, Data)

plot(pPCA$S[analysisTree$tip.label,2], pPCA$S[analysisTree$tip.label,3], pch=CladePCH[Clades], bg=CladeCols[Clades])
for (i in 1:length(CladeCols))	{
	Subset <- which(Clades == names(CladeCols)[i])
	Hull <- chull(pPCA$S[Subset,2], pPCA$S[Subset,3])
	lines(pPCA$S[c(Hull, Hull[1]),2], pPCA$S[c(Hull, Hull[1]),3], col=CladeCols[i])
}

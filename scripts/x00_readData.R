source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\cal3_Functions.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\placeFossil.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\addScale.R")

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")

library(phytools)

Tree <- read.nexus("carch_nexus.tre")
Arms <- read.csv("carch_measures.csv", stringsAsFactors = F, row.names = 1)
Ages <- read.table("carch_ages.txt")
Tree <- bind.tip(Tree, "Deinocheirus", where=grep("Bei", Tree$tip.label))

CladeVec <- rep("other theropods", nrow(Arms))
CladeVec[which(Arms[,1]=="Allosauroidea")] <- "Allosauroidea"
CladeVec[which(Arms[,1]=="Ceratosauria")] <- "Ceratosauria"
CladeVec[which(Arms[,1]=="Tyrannosauroidea")] <- "Tyrannosauroidea"
names(CladeVec) <- rownames(Arms)

Rates <- c(0.1, 0.08, 0.1) # From BAMM of Benson tree for theropods

ScaledT <- timescaleTree(Tree, Ages, Rates[1], Rates[2], Rates[3], polies = TRUE, use="ran")
addTree <- ScaledT
Drop <- c("Deltadromaeus_agilis", "Szechuanosaurus_zigongensis_(ZDM_9011)")
addTree <- drop.tip(ScaledT, Drop)

# addTree <- try(placeFossil(addTree, age = runif(1, max=77.4, min=68.6) - 66, Name = "Deinocheirus", taxa=addTree$tip.label[grep("Beishan", addTree$tip.label)], rates=Rates))
#	addTree <- try(placeFossil(addTree, age = runif(1, max=97, min=93) - 66, Name = "Meraxes", taxa=addTree$tip.label[grep("Acrocanth", addTree$tip.label)], rates=Rates))

analysisTree <- addTree
write.tree(analysisTree, "carchTree.tre")

#pdf("tree.pdf", height=5,width=8)
#plot(analysisTree, cex=0.5)
#dev.off()

#CladeCols <- setNames(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'), unique(Arms[analysisTree$tip.label,1]))
CladeCols <- setNames(c('gray60','#1b9e77','#d95f02','#7570b3'), unique(CladeVec[analysisTree$tip.label]))
Clades <- setNames(CladeVec[analysisTree$tip.label], analysisTree$tip.label)
CladePCH <- setNames(c(1, 2, 5, 8, 15, 16, 17, 3, 0, 9), unique(Arms[analysisTree$tip.label,1]))

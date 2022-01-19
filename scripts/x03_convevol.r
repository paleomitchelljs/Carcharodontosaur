library(geiger)
library(phytools)
library(convevol)

source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\cal3_Functions.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\placeFossil.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\addScale.R")

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")
source("x00_readData.R")
source("x00_convevol_Fxns.R")
source("x02_ppca.R")

# reformat in data
Data <- cbind(apply(Arms[analysisTree$tip.label,2:4], 1, sum), Arms[analysisTree$tip.label,5])
Data2 <- Arms[analysisTree$tip.label,2:5]
Data3 <- pPCA$S

out_all <- list()
outp_all <- list()
Tips_all <- c("Aucasaurus_garridoi", "Carnotaurus_sastrei", "Meraxes", "Acrocanthosaurus_atokensis", analysisTree$tip.label[grep("Tarbosaur", analysisTree$tip.label)], analysisTree$tip.label[grep("Tyrannosaur", analysisTree$tip.label)], analysisTree$tip.label[grep("Daspletosaurus", analysisTree$tip.label)], analysisTree$tip.label[grep("Gorgosaur", analysisTree$tip.label)], analysisTree$tip.label[grep("Albertosaur", analysisTree$tip.label)])

out_m <- list()
outp_m <- list()
Tips_m <- c("Meraxes", analysisTree$tip.label[grep("Tarbosaur", analysisTree$tip.label)])

out_a <- list()
outp_a <- list()
Tips_a <- c("Meraxes", "Carnotaurus_sastrei")

out_ab <- list()
outp_ab <- list()
Tips_ab <- c(analysisTree$tip.label[grep("Tarbosaur", analysisTree$tip.label)], "Carnotaurus_sastrei")

for (i in 1:100)	{
	source("x00_readData.R")

	out_all[[i]] <- convSig(analysisTree, Data, Tips_all)
	outp_all[[i]] <- convSig(analysisTree, Data3, Tips_all)

	out_m[[i]] <- convSig(analysisTree, Data, Tips_m)
	outp_m[[i]] <- convSig(analysisTree, Data3, Tips_m)

	out_a[[i]] <- convSig(analysisTree, Data, Tips_a)
	outp_a[[i]] <- convSig(analysisTree, Data3, Tips_a)

	out_ab[[i]] <- convSig(analysisTree, Data, Tips_ab)
	outp_ab[[i]] <- convSig(analysisTree, Data3, Tips_ab)
	cat("Done with run ", i, "\n")
}

CCMatrix <- list(out_all, out_m, out_a, out_ab, outp_all, outp_m, outp_a, outp_ab)
saveRDS(CCMatrix, file = "staytonCC_output.rds")

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")

Dat <- read.csv("theropod_skull.csv")

Forelimb <- Dat$Humeral.Length + as.numeric(Dat$Radial.Length) + as.numeric(Dat$MCI.length)

Ratio <- Forelimb / Dat$Femoral.length

plot(Dat$Interquadrate.width, Ratio, log="x")
text(Dat$Interquadrate.width, Ratio, Dat[,1])

KeepRows <- which(!is.na(apply(Dat[,3:ncol(Dat)],1,sum)))

Dat <- Dat[KeepRows,]
Forelimb <- Forelimb[KeepRows]
Drop <- which(Dat[,1] == "Szechuanosaurus_zigongensis_(ZDM_9011)")
Dat <- Dat[-Drop,]
Forelimb <- Forelimb[-Drop]

z <- glm(Forelimb ~ Dat$Femoral.length * Dat$Interquadrate.width)

centerscale <- function(x, LOG=T)	{
	if (LOG)	{
		x <- log(x)
	}
	x <- ( x - mean(na.omit(x)) ) / sd(na.omit(x))
	return(x)
}
Data <- data.frame(species = Dat[,1], femur = centerscale(Dat$Femoral.length), iqw = centerscale(Dat$Interquadrate.width), arm = centerscale(Forelimb))
Dat2 <- data.frame(femur = centerscale(Dat$Femoral.length), iqw = centerscale(Dat$Interquadrate.width), arm = centerscale(Forelimb))
rownames(Dat2) <- Data$species

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")
library(phytools)
library(brms)
tree <- read.tree("carchTree.tre")
DropTips <- setdiff(tree$tip.label, Data$species)
tree2 <- drop.tip(tree, DropTips)

# estimate lambda and adjust VCV matrix accordingly
A <- ape::vcv.phylo(tree2)

Mod <- bf(arm ~ femur + iqw + (1 | gr(species, cov = A)))
prior <- get_prior( Mod, data = Data, data2 = list(A = A) )

prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 1)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 1)"
prior$prior[which(prior$group == "species")] <- "normal(0, 0.01)"

#prior$prior[which(prior$coef == "femur")] <- "student_t(3, 0, 2)"
#prior$prior[which(prior$coef == "iqw")] <- "student_t(3, 0, 2)"
#prior$prior[which(prior$coef == "femur:iqw")] <- "student_t(3, 0, 2)"

ITER <- 2e4

Fit <- brm(Mod, 
		data=Data, 
		data2 = list(A=A), 
		prior = prior,
		cores = 4,
		warmup = ITER / 2,
		save_pars = save_pars(all = TRUE),
		iter = ITER,
		control = list(adapt_delta = 0.99, max_treedepth = 12)
		)


ME <- conditional_effects(Fit)

library(ggplot2)
plot(ME, points = TRUE, theme = theme_bw())

pred.1 <- predict(Fit)
plot(Data$arm[KeepRows], pred.1[,1])
plot(Data$femur[KeepRows], pred.1[,1]) 



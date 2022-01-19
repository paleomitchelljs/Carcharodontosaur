setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")
library(RRphylo)
library(phytools)
library(vioplot)

source("x00_readData.R")
source("x04_convevol_Plots.R")
source("x06_ridgeRegressionPlots.R")


pdf("threeplot_format.pdf", height=7/3, width=7, family="sans")

########## COMPOSITE PLOT
par(mar=c(2.5,3,1.25,0.75), las = 1, mfrow=c(1,3), mgp=c(1.5, 0.15, 0), tck=-0.01, cex.axis=0.67)

# PANEL A
# Stayton convergence metrics boxplots
Cmetric <- 2
boxplot(allTax[Cmetric,], merTax[Cmetric,], abelTax[Cmetric,], boxwex=0.25, col='white', names=c("","",""), ylim=c(0, 0.2), ylab="")
mtext("all", side=1, line=0.25, at=1, cex=0.5)
mtext("Meraxes / \nTarbosaurus", side=1, line=.75, at=2, font=3, cex=0.5)
mtext("Meraxes / \nCarnotaurus", side=1, line=.75, at=3, font=3, cex=0.5)
mtext(las=0, "p-val for C2", side=2, line=1.5, cex=0.75)
abline(h=0.05, lty=2, lwd=2)

mtext("A", side=3, at=0.3)

# PANEL B
# phenogram showing FEMUR on X and RATIO on Y

Xrange <- sort(nodeHeights(analysisTree)[,1], decreasing = TRUE)
Choice <- 2
Yvals <- Data[,1] / Data[,2]
ScaleBottom <- -0.1
ScaleTop <- ScaleBottom + 0.01

MapPh <- make.simmap(analysisTree, Clades, model="ER", nsim = 1)
phylomorphospace(MapPh, cbind(Data[,2], Ratio), colors = CladeCols, label = "off", xlab="femur length (mm)", ylab="log forelimb / femur length", node.size=c(0.01, 0.25))
text(Data[grep("Auca", rownames(Data)),2], Ratio[grep("Auca", names(Ratio))], "Abel.", pos=2, col=CladeCols["Ceratosauria"])
text(Data[grep("Acro", rownames(Data)),2], Ratio[grep("Acro", names(Ratio))], "Carch.", pos=3, col=CladeCols["Allosauroidea"])
text(1225, -0.35, "Tyranno.", pos=1, col=CladeCols["Tyrannosauroidea"])
points(Data[grep("Meraxes", rownames(Data)),2], Ratio[grep("Meraxes", names(Ratio))], col=CladeCols["Allosauroidea"], pch=8, cex=1.25)

#legend("top", text.col=CladeCols[1:5], legend=names(CladeCols)[1:5], bty="n")
#legend("topright", text.col=CladeCols[6:10], legend=names(CladeCols)[6:10], bty="n")
mtext("B", side=3, at=-100)

# PANEL C
# Ratio contMap from RRPHYLO
plot(Map, fsize=c(0,1), sig=2, legend=FALSE, ftype="off")
add.color.bar(90,Map$cols,title="Arm/Femur", lims=Map$lims,digits=2,prompt=FALSE,x=80, y=2,lwd=2.5,fsize=1,subtitle="")
mtext("C", side=3, at=0, line=-1.4)

dev.off()

#########
#########
#########
#########
addSil <- function (img, x = NULL, y = NULL, ysize = 1, xsize = 1, alpha = 1, color = NULL) {                                                                         
	img <- rphylopic:::recolor_phylopic(img, alpha, color)
	graphics::rasterImage(img, x - xsize/2, y - ysize/2, x + xsize/2, y + ysize/2, interpolate = TRUE)                         
} 

pdf("threeplot_sil.pdf", height=7/3, width=7, family="sans")

########## COMPOSITE PLOT
par(mar=c(2.5,3,1.25,0.75), las = 1, mfrow=c(1,3), mgp=c(1.5, 0.15, 0), tck=-0.01, cex.axis=0.67)

# PANEL A
# Stayton convergence metrics boxplots
Cmetric <- 2
boxplot(allTax[Cmetric,], merTax[Cmetric,], abelTax[Cmetric,], boxwex=0.25, col='white', names=c("","",""), ylim=c(0, 0.2), ylab="")
mtext("all", side=1, line=0.25, at=1, cex=0.5)
mtext("Meraxes / \nTarbosaurus", side=1, line=.75, at=2, font=3, cex=0.5)
mtext("Meraxes / \nCarnotaurus", side=1, line=.75, at=3, font=3, cex=0.5)
mtext(las=0, "p-val for C2", side=2, line=1.5, cex=0.75)
abline(h=0.05, lty=2, lwd=2)

mtext("A", side=3, at=0.3)

# PANEL B
# phenogram showing FEMUR on X and RATIO on Y

Xrange <- sort(nodeHeights(analysisTree)[,1], decreasing = TRUE)
Choice <- 2
Yvals <- Data[,1] / Data[,2]
ScaleBottom <- -0.1
ScaleTop <- ScaleBottom + 0.01

MapPh <- make.simmap(analysisTree, Clades, model="ER", nsim = 1)
phylomorphospace(MapPh, cbind(Data[,2], Ratio), colors = CladeCols, label = "off", xlab="femur length (mm)", ylab="log forelimb / femur length", node.size=c(0.01, 0.25))
text(Data[grep("Auca", rownames(Data)),2], Ratio[grep("Auca", names(Ratio))], "Abel.", pos=2, col=CladeCols["Ceratosauria"])
text(Data[grep("Acro", rownames(Data)),2], Ratio[grep("Acro", names(Ratio))], "Carch.", pos=3, col=CladeCols["Allosauroidea"])
text(1225, -0.35, "Tyranno.", pos=1, col=CladeCols["Tyrannosauroidea"])
points(Data[grep("Meraxes", rownames(Data)),2], Ratio[grep("Meraxes", names(Ratio))], col=CladeCols["Allosauroidea"], pch=8, cex=1.25)

#legend("top", text.col=CladeCols[1:5], legend=names(CladeCols)[1:5], bty="n")
#legend("topright", text.col=CladeCols[6:10], legend=names(CladeCols)[6:10], bty="n")
mtext("B", side=3, at=-100)

# PANEL C
# Ratio contMap from RRPHYLO
plot(Map, fsize=c(1,1), sig=2, legend=FALSE, fcolor="white", ftype="off", xlim=c(0, 230))
add.color.bar(110,Map$cols,title="Arm/Femur", lims=Map$lims,digits=2,prompt=FALSE,x=80, y=2,lwd=2.5,fsize=1,subtitle="")
mtext("C", side=3, at=0, line=-1.15)

# Carno
carn <- image_data("9cd6e394-ea75-4d75-bfb6-8ba4bc855aa3", size = 512)[[1]]
addSil(carn, x=190, y=10, ysize=4, xsize=28, alpha=1, color="black")

# Carchy
giga <- image_data("a7d7618f-f0db-47b0-b5f9-75176e7a9af5", size = 512)[[1]]
addSil(giga, 175, 15, 3.5, 32, alpha=1)

# Tyranno
tyranno <- image_data("812299bb-c476-4b9a-a918-6119f96ff1e3", size = 512)[[1]]
addSil(tyranno, 202, 23, ysize=3.7, xsize=34, alpha=1)

# Mono
mono <- image_data("e175b06d-60c1-4d0b-8596-fd2543738d19", size = 512)[[1]]
addSil(mono, 192, 41, ysize=4, xsize=22, alpha=1)

# Deino
deino <- image_data("f73a57a3-2489-44fe-9b26-d213d20ed8d0", size = 512)[[1]]
addSil(deino, 200, 34, ysize=4, xsize=28, alpha=1)

# graptor
graptor <- image_data("7495ab90-0a76-49f4-b71e-677d7bb8f7e0", size = 512)[[1]]
addSil(graptor, 197, 47, ysize=4, xsize=24, alpha=1)

# vraptor
vraptor <- image_data("224a1727-55f6-4dc0-9670-5a1d34b227cf", size = 512)[[1]]
addSil(vraptor, 195, 59, ysize=4, xsize=32, alpha=1)

dev.off()

#############################################
#############################################
#############################################
#############################################
#############################################
#############################################
pdf("threeplot_names.pdf", height=7/3, width=7, family="sans")

########## COMPOSITE PLOT
par(mar=c(2.5,3,1.25,0.75), las = 1, mfrow=c(1,3), mgp=c(1.5, 0.15, 0), tck=-0.01, cex.axis=0.67)

# PANEL A
# Stayton convergence metrics boxplots
Cmetric <- 2
boxplot(allTax[Cmetric,], merTax[Cmetric,], abelTax[Cmetric,], boxwex=0.25, col='white', names=c("","",""), ylim=c(0, 0.2), ylab="")
mtext("all", side=1, line=0.25, at=1, cex=0.5)
mtext("Meraxes / \nTarbosaurus", side=1, line=.75, at=2, font=3, cex=0.5)
mtext("Meraxes / \nCarnotaurus", side=1, line=.75, at=3, font=3, cex=0.5)
mtext(las=0, "p-val for C2", side=2, line=1.5, cex=0.75)
abline(h=0.05, lty=2, lwd=2)

mtext("A", side=3, at=0.3)

# PANEL B
# phenogram showing FEMUR on X and RATIO on Y

Xrange <- sort(nodeHeights(analysisTree)[,1], decreasing = TRUE)
Choice <- 2
Yvals <- Data[,1] / Data[,2]
ScaleBottom <- -0.1
ScaleTop <- ScaleBottom + 0.01

MapPh <- make.simmap(analysisTree, Clades, model="ER", nsim = 1)
phylomorphospace(MapPh, cbind(Data[,2], Ratio), colors = CladeCols, label = "off", xlab="femur length (mm)", ylab="log forelimb / femur length", node.size=c(0.01, 0.25))
text(Data[grep("Auca", rownames(Data)),2], Ratio[grep("Auca", names(Ratio))], "Abel.", pos=2, col=CladeCols["Ceratosauria"])
text(Data[grep("Acro", rownames(Data)),2], Ratio[grep("Acro", names(Ratio))], "Carch.", pos=3, col=CladeCols["Allosauroidea"])
text(1225, -0.35, "Tyranno.", pos=1, col=CladeCols["Tyrannosauroidea"])
points(Data[grep("Meraxes", rownames(Data)),2], Ratio[grep("Meraxes", names(Ratio))], col=CladeCols["Allosauroidea"], pch=8, cex=1.25)

#legend("top", text.col=CladeCols[1:5], legend=names(CladeCols)[1:5], bty="n")
#legend("topright", text.col=CladeCols[6:10], legend=names(CladeCols)[6:10], bty="n")
mtext("B", side=3, at=-100)

# PANEL C
# Ratio contMap from RRPHYLO
plot(Map, fsize=c(0.25,1), sig=2, legend=FALSE)
add.color.bar(90,Map$cols,title="Arm/Femur", lims=Map$lims,digits=2,prompt=FALSE,x=80, y=2,lwd=2,fsize=1,subtitle="")
mtext("C", side=3, at=0, line=-1.15)

dev.off()



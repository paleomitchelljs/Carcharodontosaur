library(phytools)

source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\cal3_Functions.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\placeFossil.R")
source("C:\\Users\\jonsm\\Dropbox\\Research\\BAMMtools_Functions\\addScale.R")

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\CarcharoCampanas")
source("x00_readData.R")

########## 
# ratio
Data <- cbind(apply(Arms[analysisTree$tip.label,2:4], 1, sum), Arms[analysisTree$tip.label,5])
Map <- make.simmap(analysisTree, Clades, model="ER", nsim = 1)

Xrange <- sort(nodeHeights(analysisTree)[,1], decreasing = TRUE)
Choice <- 2
Yvals <- Data[,1] / Data[,2]
ScaleBottom <- -0.1
ScaleTop <- ScaleBottom + 0.01

MaxAge <- max(nodeHeights(analysisTree))

pdf("ratio_labels.pdf", height=5, width=6.5)
par(xaxt="n", mar=c(4,4,1,1), las = 1)
locs <- phenogram(Map, Yvals, spread.labels=F, fsize = 0.5, xlab="", ylab="forelimb / femur", lwd=2, colors = CladeCols)

par(new=TRUE, xaxt="s", mgp=c(1.75, 0.5, 0), tck=-0.01)
plot(1, type="n", xlab="", ylab="", xlim=c(MaxAge+66, 65), ylim=c(0,1), yaxt="n", xaxt="n", xpd=NA)
addScale(Xrange + 66, YLim = c(ScaleBottom, ScaleTop))
axis(1, at=c(500, ceiling(MaxAge+66), 201, 145, 66))
mtext(side = 1, text = "age (Ma)", line = 1.5)
legend("topleft", pt.bg = CladeCols, legend=names(CladeCols), bty="n", pch=21, pt.cex=1.5)
dev.off()

MerAge <- nodeHeights(analysisTree)[which(analysisTree$edge[,2] == which(analysisTree$tip.label == "Meraxes")),2]
pdf("ratio_nolabels.pdf", height=5, width=6.5)
par(xaxt="n", mar=c(4,4,1,1), las = 1)
locs <- phenogram(Map, Yvals, spread.labels=F, fsize = 0, xlab="", ylab="forelimb / femur", lwd=2, colors = CladeCols)
arrows(x0 = MerAge, y0 = 0.25, x1 = MerAge, y1 = Yvals["Meraxes"]-0.03, lwd=2, angle = 25, length = 0.15)

par(new=TRUE, xaxt="s", mgp=c(1.75, 0.5, 0), tck=-0.01)
plot(1, type="n", xlab="", ylab="", xlim=c(MaxAge, 65), ylim=c(0,1), yaxt="n", xaxt="n", xpd=NA)
addScale(Xrange + 66, YLim = c(ScaleBottom, ScaleTop))
axis(1, at=c(500, ceiling(MaxAge+66), 201, 145, 66))
mtext(side = 1, text = "age (Ma)", line = 1.5)
legend("topleft", pt.bg = CladeCols, legend=names(CladeCols), bty="n", pch=21, pt.cex=1.5)
dev.off()

pdf("theropod_phylomorph.pdf", height=6.5, width=6.5)
par(mar=c(4,4,1,1), las = 1, mgp=c(1.75, 0.5, 0), tck=-0.01)
phylomorphospace(Map, cbind(Data[,2], Yvals), colors = CladeCols, label = "off", xlab="femur length (mm)", ylab="forelimb / femur length")
arrows(x0 = Data["Meraxes",2], y0 = 0.25, x1 = Data["Meraxes",2], y1 = Yvals["Meraxes"]-0.03, lwd=2, angle=25)
dev.off()


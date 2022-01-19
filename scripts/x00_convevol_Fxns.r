#conevol3
# find vector of ancestors
pullNodeSeq <- function(phy, t1)	{
	if (is.na(as.numeric(t1)))	{
		t1 <- which(phy$tip.label == t1)
	}
	t1 <- as.numeric(t1)
	Stop <- Ntip(phy) + 1
	Anc <- phy$edge[which(phy$edge[,2] == t1),1]
	if (Anc != Stop)		{
		Anc <- c(Anc, pullNodeSeq(phy, Anc))
	}
	return(as.character(Anc))
}
# calc distance between ancestors of two tips and their ancestors
calcCs <- function(tips, ancList, allDists, phy, VERBOSE=FALSE)	{
	if (VERBOSE == TRUE)		{
		cat("Analzying ", tips)
	}
	# tip distances
	Dtip <- allDists[tips[1], tips[2]]
	
	# what nodes are shared between tip1 and tip2?
	overLap <- intersect(ancList[[tips[1]]], ancList[[tips[2]]])
	
	# what node is the mrca of tip1 and tip2?
	MRCA <- max(as.numeric(overLap))
	
	# what nodes are UNIque to tip 1 and what are UNIque to tip 2? E.g., what nodes lead from the MRCA to each tip?
	unit1 <- c(setdiff(c(tips[1], ancList[[tips[1]]]), overLap), MRCA)
	unit2 <- c(setdiff(c(tips[2], ancList[[tips[2]]]), overLap), MRCA)
	
	# compare all post-MRCA nodes along each lineage to one another (including the tips)
	ancPairs <- sapply(unit1, function(t1) sapply(unit2, function(t2) allDists[t1,t2]))
	
	# sum changes along each lineage including the MRCA
	lin1 <- c(unit1, as.character(MRCA))
	Lin1 <- sum(sapply(1:(length(lin1)-1), function(x) allDists[lin1[x],lin1[x+1]]))
	lin2 <- c(unit2, as.character(MRCA))
	Lin2 <- sum(sapply(1:(length(lin2)-1), function(x) allDists[lin2[x],lin2[x+1]]))
	
	# sum all pairwise distances between nodes of the subtree from the MRCA of tip1 & tip2
	totalClade <- as.character(getDescendants(phy, MRCA))
	totalMove <- sum(allDists[totalClade, totalClade]) / 2
	
	# Dmax as the maximum distance between ancestor pairs
	Dmax <- max(ancPairs)
	C1 <- 1 - (Dtip / Dmax)
	C2 <- Dmax - Dtip
	C3 <- C2 / (Lin1 + Lin2)
	C4 <- C2 / totalMove
	return(c(C1=C1, C2=C2, C3=C3, C4=C4))
}

# workhorse function
calcConv <- function(phy, traits, focaltaxa, VERBOSE=FALSE, Method="bm")	{
	require(phytools)
	if (class(phy) != "phylo")	{
		stop("first argument, phy, must be of class phylo!")
	}
	if (class(traits) != "matrix")	{
		Names <- names(traits)
		traits <- matrix(traits, ncol=1)
		rownames(traits) <- Names
	}
	if (nrow(traits) != Ntip(phy))	{
		stop("Number of rows does not match number of tips!")
	}
	if (length(intersect(rownames(traits), phy$tip.label)) != Ntip(phy))		{
		print("Warning: the row names in traits do not match tip labels in phy; assume that they are in the same order.")
		rownames(traits) <- phy$tip.label
	}
	traits <- traits[phy$tip.label,]
	if (VERBOSE == TRUE)		{
		print("reached traits")
	}
	# for C1: need phylogeny, trait matrix, ancestral recos
	if (Method == "bm")	{
		ACE <- apply(traits, 2, fastAnc, tree=phy)
	}
	else if (Method == "rr")	{
		require(RRphylo)
		ACE <- RRphylo()
		error("TKTKTK Not implemented yet")
	}
	# number the tips
	tipNums <- as.character(1:Ntip(phy))
	names(tipNums) <- phy$tip.label
	focaltaxa <- tipNums[focaltaxa]
	rownames(traits) <- tipNums[rownames(traits)]

	# matrix with node and tip vectors
	allVals <- rbind(traits, ACE)

	# get distances between tips and nodes
	allDists <- as.matrix(dist(allVals, method="euclidean"))

	# What node is ancestral to each tip?
	Ancs <- lapply(tipNums, pullNodeSeq, phy=phy)
	names(Ancs) <- tipNums

	# unique combinations of focal taxa
	Combinations <- combn(focaltaxa, 2)
	if (VERBOSE == TRUE)		{
		print("starting combinations...")
	}
	# calculate the C values 1-4
	Cmat <- apply(Combinations, 2, function(x) calcCs(x, ancList=Ancs, allDists=allDists, phy=phy, VERBOSE=VERBOSE))

	out <- apply(Cmat,1,mean)
	return(out)
}

convSig <- function(phy, traits, focaltaxa, nsim=1e2)	{
	data <- calcConv(phy, traits, focaltaxa)
	
	# run simulations
	phylMat <- vcv.phylo(phy)
	phylMat2 <- phyl.vcv(traits, phylMat, 0)
	simDat <- sim.char(phy, phylMat2$R, nsim, model="BM", root=0)
	simOut <- apply(simDat, 3, calcConv, phy=phy, focaltaxa=focaltaxa)
	pvals <- sapply(1:4, function(x) length(which(simOut[x,] >= data[x]))) / nsim
	out <- cbind(data, pvals)
	return(out)
}

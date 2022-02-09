#' @title Extract Fossil Record From Tree

#' @description
#' Thus function collects the fossil record (ages) of a particular clade from the result simmulations in the package [FossilSim]

#' @param x the node number that defines the crown-clade
#' @param phy a phylogenetic tree of class 'phylo'
#' @param fossils an object of class 'fossils' produced by [FossiSim]

#' @details
#' Coming soon.

#' @return a vector of fossil ages

#' @import phangorn

#' @export

fossil.record <- function(x, phy, fossils) {

	# Phangorn function Descendants identifies the descendant nodes and tips
	
	des <- phangorn::Descendants(phy, node=x, type = "all")

	# Then use the descendants indices to retrieve the fossil ages from the fossil object
	fages <- numeric()

	for(i in 1:length(des)) {
	fages <- c(fages, fossils$hmin[which(fossils$edge == des[i])])
	}
	
	# Get tip names to faciliate the generation of the clade constraint in MrBayes

	tnames <- phy$tip.label[Descendants(phy, node=x, type = "tips")[[1]]]

	RES <- list(mrca.node = x, tip.names=tnames, fossil.ages=fages, n.foss=length(fages))
	
	return(RES)
}


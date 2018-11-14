library(CellNOptR)

# test readBND ------------------------------------------------------------------

required_pcks = list("plyr","dplyr","tidyr","readr")

if(!all(unlist(lapply(required_pcks,requireNamespace)))){
	print("the following packages need to be installed to use readBND:")	
	print(unlist(required_pcks))
	print("Please, install the packages manually for this feature.")

}else{
	#download.file("https://maboss.curie.fr/pub/example.bnd",destfile = "./tests/example.bnd")
	#model = readBND("./tests/example.bnd")
	model = readBND("https://maboss.curie.fr/pub/example.bnd")
	
	# basic checks for being a CellNoptR model:
	stopifnot(is.list(model))
	stopifnot(length(model)==4)
	stopifnot(all(names(model) %in% c("reacID","namesSpecies","interMat","notMat")))
	stopifnot(length(model$reacID)==ncol(model$interMat))
	stopifnot(length(model$reacID)==ncol(model$notMat))
	stopifnot(length(model$namesSpecies)==nrow(model$interMat))
	stopifnot(length(model$namesSpecies)==nrow(model$notMat))
	
	#plotModel(model)
	
	
	
}



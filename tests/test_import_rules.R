library(CellNOptR)

# test readBND ------------------------------------------------------------------

?readBND
download.file("https://maboss.curie.fr/pub/example.bnd",destfile = "./tests/example.bnd")

model = readBND("./tests/example.bnd")

# basic checks for being a CellNoptR model:
stopifnot(is.list(model))
stopifnot(length(model)==4)
stopifnot(all(names(model) %in% c("reacID","namesSpecies","interMat","notMat")))
stopifnot(length(model$reacID)==ncol(model$interMat))
stopifnot(length(model$reacID)==ncol(model$notMat))
stopifnot(length(model$namesSpecies)==nrow(model$interMat))
stopifnot(length(model$namesSpecies)==nrow(model$notMat))

#plotModel(model)

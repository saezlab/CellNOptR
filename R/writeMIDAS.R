writeMIDAS <- function(CNOlist, filename)
{

    cnolist = CNOlist

    namesCues = cnolist$namesCues
    namesStimuli = cnolist$namesStimuli
    namesInhibitors = cnolist$namesInhibitors
    namesSignals = cnolist$namesSignals
    timeSignals = cnolist$timeSignals

    nCues = length(namesCues)
    nStimuli = length(namesStimuli)
    nInhibitors = length(namesInhibitors)
    nSignals = length(namesSignals)
    nTimes = length(timeSignals)

    # First column is cellline
    nCols = 1 + nStimuli + nInhibitors + nSignals * 2 
    nConds = dim(cnolist$valueSignals[[1]])[1]
    nRows = nConds * nTimes

    # build the header
    colNames = c("TR:mock:CellLine")
    for (x in namesStimuli){
        colNames = cbind(colNames, paste("TR:", x, sep=""))
    }
    for (x in namesInhibitors){
        colNames = cbind(colNames, paste("TR:", x, "i", sep=""))
    }
    for (x in namesSignals){
        colNames = cbind(colNames, paste("DA:", x, sep=""))
    }
    for (x in namesSignals){
        colNames = cbind(colNames, paste("DV:", x, sep=""))
    }



    data = matrix(NA, nRows, nCols)

    # cellline is made of ones for now
    data[,1] = 1

    # fill the experiments first
    for (time in 1:nTimes){

        i1 = 1 + (time - 1) * nConds
        i2 = i1 + nConds - 1

        j1 = 2; j2 = j1 + nStimuli - 1
        data[i1:i2, j1:j2] = cnolist$valueStimuli

        j1 = j2 + 1; j2 = j1 + nInhibitors - 1
        data[i1:i2, j1:j2] = cnolist$valueInhibitors

        j1 = j2 + 1; j2 = j1 + nSignals - 1
        data[i1:i2, j1:j2] = rep(timeSignals[time], nSignals)

        j1 = j2 + 1; j2 = j1 + nSignals - 1
        data[i1:i2,j1:j2] = cnolist$valueSignals[[time]] 
    }


    # and save it into a file
    if (file.exists(filename)==FALSE){
        write.table(data, file=filename,
            row.names=FALSE, col.names=colNames,quote=FALSE,sep=",")
    }
    else{
       stop(paste("File ", filename, "already exists.",  sep=""))
    }





}



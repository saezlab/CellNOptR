#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$
makeCNOlist<-function(dataset,subfield, verbose=TRUE){

    #check that all the needed elements are present
    if(!is.list(dataset)){
        stop("The input to this function should be a list with elements 'dataMatrix', 'TRcol','DAcol',and 'DVcol'")
        }

    if(length(dataset) != 4){
        stop("The input to this function should be a list with elements 'dataMatrix', 'TRcol','DAcol',and 'DVcol'")
        }

    if(sum(c("dataMatrix","TRcol","DAcol","DVcol") %in% names(dataset))!=4){
        stop("The input to this function should be a list with elements 'dataMatrix', 'TRcol','DAcol',and 'DVcol'")
        }

    # first, we summarise the replicates in the dataMatrix: replicates are rows in the dataset$dataMatrix
    # that have the exact same values in all TR: and DA: columns
    duplCond<-as.matrix(dataset$dataMatrix[,c(dataset$TRcol,dataset$DAcol)])
    duplRows<-which(duplicated(duplCond) == TRUE)
    if (verbose == TRUE){
        if (length(duplRows)>0){
            print("Cleaning duplicated rows")
        }
    }

    # creates a variance matrix
    variances = dataset$dataMatrix * 0


    while(length(duplRows) != 0){
        # the all(x == ) is buggy in the case of NA hence the compareNA function.
        #dupIndex<-apply(duplCond,MARGIN=1,function(x) all(x == duplCond[duplRows[1],]))
        dupIndex<-apply(duplCond,MARGIN=1,function(x)all(compareNA(x,duplCond[duplRows[1],])))

        dupIndex<-which(dupIndex == TRUE)
        dupMatrix<-dataset$dataMatrix[dupIndex,]
        #compute the new row as the average across duplicate rows
        newRow<-colMeans(dupMatrix, na.rm=TRUE)
        # variance for these rows
        newVariance = apply(dupMatrix, MARGIN=2, FUN=var, na.rm=T)

        #replace the first duplicated row by the summarised one
        dataset$dataMatrix[dupIndex[1],]<-newRow

        # same for the variance
        variances[dupIndex[1],]<-newVariance

        #remove the other summarised rows
        dataset$dataMatrix<-dataset$dataMatrix[-dupIndex[2:length(dupIndex)],]
        variances<-variances[-dupIndex[2:length(dupIndex)],]

        duplCond<-as.matrix(dataset$dataMatrix[,c(dataset$TRcol,dataset$DAcol)])
        duplRows<-which(duplicated(duplCond) == TRUE)
    }

    if (any(is.nan(as.matrix(dataset$dataMatrix)))){
        dataset$dataMatrix[is.nan(as.matrix(dataset$dataMatrix))] <- NA
    }

    #now extract the names of the cues, and the inhibitors/stimuli
    namesCues<-colnames(dataset$dataMatrix)[dataset$TRcol]

    if(subfield == TRUE){
        namesCues<-sub(pattern="(TR:)",x=namesCues, replacement="",perl=TRUE)
        tagInhib<-grep(pattern=":Inhibitor", x=namesCues)
        # tagInhibL is a logical version of the previous statement. 
        # must be set now before namesCues is changed. See JIRA bug 27
        tagInhibL<-grepl(pattern=":Inhibitor", x=namesCues)

        # remove the trailing :Inhibitors and :Stimuli
        namesCues<-sub(pattern="(:\\w*$)",x=namesCues, replacement="",perl=TRUE)

        # remove trailing i
        namesCues[tagInhib]<-sub(pattern="(i$)", x=namesCues[tagInhib], replacement="", perl=TRUE)

        # if no inhibitors, grep returns integer(0), so we now need to use grepl
        # (logical version of grep)
        namesStimuli<-namesCues[tagInhibL==FALSE]
        namesInhibitors<-namesCues[tagInhibL==TRUE]
        }

    if(subfield == FALSE){
        namesCues<-sub(pattern="(TR:)",x=namesCues, replacement="",perl=TRUE)
        tagInhib<-grep(pattern="(i$)", x=namesCues, perl=TRUE,ignore.case=FALSE)
        # tagInhibL must be set now before namesCues is changed See JIRA bug 27
        tagInhibL<-grepl(pattern="(i$)", x=namesCues, perl=TRUE,ignore.case=FALSE)
        namesCues[tagInhib]<-sub(pattern="(i$)", x=namesCues[tagInhib], replacement="", perl=TRUE)
        # if no inhibitors, grep returns integer(0), so we now need to use grepl
        # (logical version of grep)
        namesStimuli<-namesCues[tagInhibL==FALSE]
        namesInhibitors<-namesCues[tagInhibL==TRUE]

    }

    if(sum("NOCYTO" %in% namesCues) != 0){
        namesCues<-namesCues[-grep(pattern="NOCYTO", namesCues)]
        namesStimuli<-namesStimuli[-grep(pattern="NOCYTO", namesStimuli)]
        }

    if(sum("NOINHIB" %in% namesCues) != 0){
        namesCues<-namesCues[-grep(pattern="NOINHIB", namesCues)]
        namesStimuli<-namesStimuli[-grep(pattern="NOINHIB", namesStimuli)]
        }

    if(sum("NO-INHIB" %in% namesCues) != 0){
        stop("Found a column with NO-INHIB tag. MIDAS files must use NOINHIB instead. Fix your MIDAS file please")
    }
    if(sum("NO-LIG" %in% namesCues) != 0){
        stop("Found a column with NO-LIG tag. MIDAS files do not accept NO-LIG. use NOINHIB or NOCYTO instead. Fix your MIDAS file please")
    }
    if(sum("NOLIG" %in% namesCues) != 0){
        stop("Found a column with NO-LIG tag. MIDAS files do not accept NO-LIG. use NOINHIB or NOCYTO instead. Fix your MIDAS file please")
    }
    if(sum("NO-CYTO" %in% namesCues) != 0){
        stop("Found a column with NO-CYTO tag. MIDAS file must use NOCYTO instead. Fix your MIDAS file please")
    }

    #now extract the names of the signals
    namesSignals<-colnames(dataset$dataMatrix)[dataset$DAcol]
    namesSignals<-sub(
        pattern="(DA:p-)",
        x=namesSignals,
        replacement="",
        perl=TRUE)
    namesSignals<-sub(
        pattern="(DA:)",
        x=namesSignals,
        replacement="",
        perl=TRUE)

    #now extract the time signals
    times<-as.factor(as.vector(as.character(as.matrix(dataset$dataMatrix[,dataset$DAcol]))))
    timeSignals<-sort(as.double(levels(times)))

    #Build the valueCues matrix (i.e. a matrix with nrows=nrows in dataMatrix and ncol=number of cues,
    #filled with 0/1 if the particular cue is present or not)

#1.I create a matrix that is a subset of the data, and only contains the TR columns
#(the cellLine TR column was removed previously)

#2.I remove the columns with NOCYTO or NOINHIB (if they exist), they don't bring any info

    if(length(grep(pattern="NOCYTO",colnames(dataset$dataMatrix)[dataset$TRcol])) != 0){
        nocyto<-grep(pattern="NOCYTO",colnames(dataset$dataMatrix)[dataset$TRcol])
        TRcol<-dataset$TRcol[-nocyto]
        cues<-dataset$dataMatrix[,TRcol]

        }else{
            # use as.matrix and then set colnames to cope for particular case of
            # only one cue
            cues<-as.matrix(dataset$dataMatrix[,dataset$TRcol])
            colnames(cues) = colnames(dataset$dataMatrix)[dataset$TRcol]
            }

    if(length(grep(pattern="NOINHIB",colnames(cues))) != 0){
        noinhib<-grep(pattern="NOINHIB",colnames(cues))
        cues<-cues[,-noinhib]
        }

#3. The cues sub-data frame now contains 1 if the cue is present and 0/NA otherwise,
#so I just need to transform the data frame into a numerical matrix
#and replace the NA in there by zeros
    cues<-as.matrix(cues,nrow=dim(cues)[1],ncol=dim(cues)[2],byrow=TRUE)
    cues[is.na(cues)]<-0

#Build the valueSignals matrices. I am going to build one big matrix
#that includes all the time points, and then I will split it into one matrix for each time point
#And then I will arrange the valueCues matrix accordingly
    valueSignals<-as.matrix(dataset$dataMatrix[,dataset$DVcol])
    valueVariance<-as.matrix(dataset$dataMatrix[,dataset$DVcol])

#This bit will create an index that contains all rows with timept1, 2, 3,...

#1.First I create a matrix that holds the time information for each row
    times<-as.matrix(dataset$dataMatrix[,dataset$DAcol])

#2. Now I check that all the columns are tha same, i.e. that each row
#will contain data on the same time point

    if (length(dataset$DAcol)>1){
        check<-rep(FALSE,(length(dataset$DAcol)-1))
        for(i in 1:length(check)){
            check[i]<-all.equal(times[,i],times[,(i+1)])
            }
        if(sum(check) != length(check))    {
           warning("Each row of your data file should contain measurements at the same time point.
               The times for the first DA column will be considered as the times for all measurements")
           }
    }
#3.Now I will only use the first column of times

#First, I create a vector timeRows that contains the indexes of the rows that contain data
#about each time point (in increasing order of time), and the vector whereTimes that contain`
#the info about how many rows I have for each time (which will allow me to extract the right
#entries from timesRows)

    times<-times[,1]
    ntimes<-length(timeSignals)
    if (ntimes <2){
        stop("Error while parsing the data. Only one time was found.")
    }
    whereTimes<-rep(0,ntimes)
    timesRows<-0  # Melody uses timesRows <- rep(0,1)
    for(i in 1:ntimes){
        timesRows<-c(timesRows,which(times == timeSignals[i]))
        whereTimes[i]<-length(which(times == timeSignals[i]))
        }

    timesRows<-timesRows[2:length(timesRows)]

#Check that we have data across all conditions for all time points except zero
    if(length(unique(whereTimes[2:length(whereTimes)])) != 1){
        warning("This program expects data across all conditions at all time points (except t=0) ")
        }
    if (verbose){
        print("Please be aware that if you only have some conditions at time zero (e.g.only inhibitor/no inhibitor), the measurements for these conditions will be copied across matching conditions at t=0")
    }
#Do the t=0 matrix, and produce a new cues matrix, that does not contain duplicates
#(1 row per condition and different matrices will be build for the different times)
    valueSignals<-list(matrix(data=0,nrow=whereTimes[2],ncol=length(dataset$DVcol)))
    valueVariance<-list(matrix(data=0,nrow=whereTimes[2],ncol=length(dataset$DVcol)))

#This vector tells me which columns of the cues matrix I should pay attention to when
#copying data across for time=0

    # bug report 31 and
    if (dim(cues)[2] >1){
        # bug report 44
        if (whereTimes[1] == 1){
            zerosCond <- 0
        }
        else{
            zerosCond<-apply(cues[timesRows[1:whereTimes[1]],],1,function(x) which(x > 0))
        }
    }
    else{
        warning("unusual case with 1 dimension in cues")
        zerosCond<-which(cues[timesRows[1:whereTimes[1]],] > 0)
    }
    zerosCond<-unique(unlist(zerosCond))
    count=1
    newcues<-matrix(data=0,nrow=whereTimes[2],ncol=dim(cues)[2])

    #fix bug report 38 to be able to have mixed times in a MIDAS file
    #for(i in timesRows[1]:timesRows[whereTimes[1]]){
    for(i in timesRows[1:whereTimes[1]]){
        # 15.12.2017: is this a bug? "i" is already going through some element of timesRow variable. 
    	# is this "timesRows[i]" make sense?! or we can just use i instead
    	#present<-zerosCond[which(cues[timesRows[i],zerosCond] > 0)]
    	present<-zerosCond[which(cues[i ,zerosCond] > 0)]
    	 
        if(length(present) == 0){
            for(n in timesRows[(whereTimes[1]+1):(whereTimes[1]+whereTimes[2])]){
                    if(sum(cues[n,zerosCond]) == 0){

                        valueSignals[[1]][count,]<-as.numeric(dataset$dataMatrix[i,dataset$DVcol])
                        valueVariance[[1]][count,]<-as.numeric(variances[i,dataset$DVcol])
                        newcues[count,]<-cues[n,]
                        count=count+1
                        }
                    }
        }else{
        	
            for(n in timesRows[(whereTimes[1]+1):(whereTimes[1]+whereTimes[2])]){
                if(length(zerosCond[which(cues[n,zerosCond] > 0)]) == length(present)){
                    if(all(zerosCond[which(cues[n,zerosCond] > 0)] == present) &&   # same cues are there
                        length(which(cues[n,zerosCond] > 0)) != 0 &&  # there is at least one non-zero cue
                       # same principle as above: use i instead of timesRows[i]
                       #all(cues[timesRows[i],zerosCond] == cues[n,zerosCond]) # cues have the same level
                       all(cues[i,zerosCond] == cues[n,zerosCond]) # cues have the same level
                       ){
                    	valueSignals[[1]][count,]<-as.numeric(dataset$dataMatrix[i,dataset$DVcol])
                        valueVariance[[1]][count,]<-as.numeric(variances[i,dataset$DVcol])
                        newcues[count,]<-cues[n,]
                        count=count+1
                    }
                }
            }
        }
    }

    #Now build the matrices for the other time points
    for(i in 2:length(timeSignals)){
        valuesTi<-matrix(data=0,nrow=whereTimes[2],ncol=length(dataset$DVcol))
        valuesVarianceTi<-matrix(data=0,nrow=whereTimes[2],ncol=length(dataset$DVcol))
        for(n in 1:dim(newcues)[1]){
            rowsMatchCues<-apply(cues,1,function(x) all(x == newcues[n,]))
            rowsmatchTime<-times == timeSignals[i]
            rowsMatch<-which((rowsMatchCues + rowsmatchTime) == 2)
            valuesTi[n,]<-as.numeric(dataset$dataMatrix[rowsMatch,dataset$DVcol])
            valuesVarianceTi[n,]<-as.numeric(variances[rowsMatch,dataset$DVcol])
            }

        valueSignals[[i]]<-valuesTi
        valueVariance[[i]]<-valuesVarianceTi

        }

#Build the valueInhibitors and valueStimuli  matrices, which are a subset of the cues one

    if(subfield == TRUE){
        valueInhibitors<-newcues[,grep(
            pattern="Inhibitor",x=colnames(cues),ignore.case=TRUE)]
        valueStimuli<-newcues[,grepl(pattern="Inhibitor",x=colnames(cues),ignore.case=TRUE)==FALSE]

        }else{

            valueInhibitors<-newcues[,grep(pattern="(i$)",x=colnames(cues),ignore.case=FALSE,perl=TRUE)]
            valueStimuli <- newcues[,grepl(pattern="(i$)",x=colnames(cues),ignore.case=FALSE,perl=TRUE)==FALSE]

            }
    if(is.null(dim(valueInhibitors))){
        valueInhibitors<-matrix(valueInhibitors,nrow=dim(newcues)[1])
        }

    if(is.null(dim(valueStimuli))){
        valueStimuli<-matrix(valueStimuli,nrow=dim(newcues)[1])
        }


    return(list(
        namesCues=namesCues,
        namesStimuli=namesStimuli,
        namesInhibitors=namesInhibitors,
        namesSignals=namesSignals,
        timeSignals=timeSignals,
        valueCues=newcues,
        valueInhibitors=valueInhibitors,
        valueStimuli=valueStimuli,
        valueSignals=valueSignals,
        valueVariances=valueVariance
    ))

}


compareNA <- function(v1,v2) {
    # This function returns TRUE wherever elements are the same, including NA's,
    # and false everywhere else.
    same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
}

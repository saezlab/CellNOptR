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
# $Id: readMIDAS.R 3464 2013-04-05 12:26:12Z cokelaer $
readMIDAS<-function(MIDASfile, verbose=TRUE){

    # Read the data. 
    # First, let us store the length of each row by counting the column
    ncols = count.fields(MIDASfile, sep = ",")

    # check that there is at least a header and one row of data (i.e. 2 rows)
    if (length(ncols) == 1){
        stop("Your input file has a wrong format (only 1 column)")
    }
    # if at least one row of data, all of them must have the same length
    if (length(unique(ncols[2:length(ncols)]))>1){
        print(unique(ncols[2:length(ncols)]))
        stop("Your input file has a wrong format (rows have different lengths)")
    }

    # Second, let us check the consistency between row and header lengths
    maxncol <- max(ncols)
    ncol_header = ncols[1]
    if (ncol_header == maxncol){ #everything seems right, let us use read.csv 
        data<-read.csv(file=MIDASfile,header=TRUE,sep=',',fill=TRUE,as.is=TRUE,check.names=FALSE)
    }
    else{
        # too much columns but we can still proceed. This is useful when extra empty
        # commas are used. We cannot use read.csv that expect length of each row
        # to correspond to the length of the header.  
        if (ncol_header < maxncol){
            warning("Some rows have more columns than the header. Extra columns will be removed. Check and fix your MIDAS file.")

            # read the data (except header) and convert to data.frame
            data <- data.frame(read.table(MIDASfile,  skip=1,sep=",",flush=TRUE))

            # remove the extra columns as compared to the header
            data <- data[1:nrow(data), 1:ncol_header]
            # finally set the dataframe column names to be those found in the header
            # use read.csv setting nrow=0 . it reads all the data (not optimal)
            # but do not care about discrepancy between columns length.
            # Note also that the check.names is required to keep track of the special keyword :
            # (otherwise replaced by a dot .) 
            header_column_names = read.csv(MIDASfile,sep=",", nrow=0, check.names=FALSE)

            colnames(data) <- colnames(header_column_names)

            #data<-read.csv(file=MIDASfile,sep=',',fill=TRUE,as.is=TRUE,col.names=paste("V",seq_len(ncol_header), sep = ""))
        }
        else {
            # not enough data as compared to the header, so let us stop. could
            # fill with NA but there is surely a mistake in such case.
            stop("Some rows have less columns than the header. Fix your MIDAS file.")
        }
    }

    # Determine which are the informative columns (i.e. columns with useful data info and values)
    TRcol<-grep(pattern="TR:",x=colnames(data),ignore.case=FALSE)
    DAcol<-grep(pattern="DA:",x=colnames(data),ignore.case=FALSE)
    DVcol<-grep(pattern="DV:",x=colnames(data),ignore.case=FALSE)
    #data<-data[,c(TRcol,DAcol,DVcol)]

    # Print information about the data set 
    if (verbose){
        print(paste(
            "Your data set comprises ", nrow(data),
            "conditions (i.e. combinations of time point and treatment)"))
    }

    # Check that the right number of columns are present            


    # The MIDAS may contain the DA:ALL special time 
    if(length(grep("DA:ALL",colnames(data))) != 0 && length(DAcol) == 1){
        if (verbose){
            print("Your data comprises only a DA:ALL column;  all readouts are assumed to have been acquired at the same time.")
        }

        # Make the additional columns
        newDA<-matrix(rep(data[,DAcol],length(DVcol)),byrow=FALSE,ncol=length(DVcol))
        newDAnames<-colnames(data)[DVcol]
        colnames(newDA)<-sub(pattern="DV:",x=newDAnames,replacement="DA:",ignore.case=FALSE)
        data<-cbind(data,newDA)
        data<-data[,-DAcol]
    } else {
        if(length(DAcol) != length(DVcol)){
            warning("DA columns and DV columns do not match.")
        }
    }

    if (verbose){
        print(paste("Your data set comprises measurements on ", length(DVcol)," different species"))
    }

    # first, replace all CellType string into CellLine.
    CellType<-grep(pattern="(TR:\\w*:CellType)",
        x=colnames(data),ignore.case=TRUE,perl=TRUE,value=TRUE)
    for (x in CellType){
        names(data)[names(data)==x] <- sub(":CellType", ":CellLine",x)
    }

    # then look at all CellLine strings before removing them.
    CellLine<-grep(pattern="(TR:\\w*:CellLine)",
        x=colnames(data),ignore.case=TRUE,perl=TRUE,value=TRUE)


    if(length(CellLine) != 0){

        CellLine<-sub(pattern="TR:",x=CellLine,replacement="",ignore.case=FALSE)
        CellLine<-sub(pattern=":CellLine",x=CellLine,replacement="",ignore.case=TRUE)

        if (verbose){
            print(paste(
                "Your data set comprises ", (length(TRcol)-length(CellLine)),
                "stimuli/inhibitors and", length(CellLine),"cell line(s) (",CellLine,")" ))
        }

        data<-data[,-grep(pattern="(TR:\\w*:CellLine)",x=colnames(data),ignore.case=TRUE,perl=TRUE,value=FALSE)]
    }
    else{
        if (verbose){
            print(paste("Your data set comprises ", length(TRcol),"stimuli and inhibitors"))
        }
        warning("There is no cell line information. If some of your TR columns represents the cell lines, please indicate it in your file by naming them 'TR:name:CellLine' (you may use TR:name:CellType' as well.")
    }

    if (verbose){
        print("Please be aware that CNO only handles measurements on one cell line at this time.")
    }

    # data has been changed so we need to extract columns indices again
    TRcol<-grep(pattern="TR:",x=colnames(data),ignore.case=FALSE)
    DAcol<-grep(pattern="DA:",x=colnames(data),ignore.case=FALSE)
    DVcol<-grep(pattern="DV:",x=colnames(data),ignore.case=FALSE)
    data<-data[,c(TRcol,DAcol,DVcol)]


    # replace NaN character by NA. as.matrix is required to scan all columns AND rows
    # Note that on some older R version, the is.nan does not seem to work well, hence the
    # try catch (see Changelog 0.99.6
    conversion <- tryCatch({
        if (any(is.nan(as.matrix(data)))){
            data[is.nan(as.matrix(data))] <- NA
        }
        conversion = TRUE},
        error=function(e) {return(FALSE)},
        finally={
            if (verbose){
                print("Your data file contained 'NaN'. We have assumed that these were missing values and replaced them by NAs.")
            }}
    )

    if ( conversion == FALSE){
        if(any(as.matrix(data == "NaN"))){
            for(c in 1:dim(data)[2]){
                for(r in 1:dim(data)[1]){
                     if(data[r,c] == NA) {
                         data[r,c]<-NA
                     }
                }
            }
            if (verbose){
                print("Your data file contained 'NaN'. We have assumed that these were missing values and replaced them by NAs.")
            }
        }
    }



    return(list(
        dataMatrix=data,
        TRcol=grep(pattern="TR:",x=colnames(data),ignore.case=FALSE),
        DAcol=grep(pattern="DA:",x=colnames(data),ignore.case=FALSE),
        DVcol=grep(pattern="DV:",x=colnames(data),ignore.case=FALSE)))

    }


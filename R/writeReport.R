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
# $Id: writeReport.R 3155 2013-01-09 15:24:58Z cokelaer $
writeReport<-function(
	modelOriginal,
	modelOpt,
	optimResT1,
	optimResT2,
	CNOlist,
	directory,
	namesFiles=list(
		dataPlot=NA,
		evolFitT1=NA,
		evolFitT2=NA,
		simResultsT1=NA,
		simResultsT2=NA,
		scaffold=NA,
		scaffoldDot=NA,
		tscaffold=NA,
		wscaffold=NA,
        PKN=NA,
		PKNdot=NA,
		wPKN=NA,
		nPKN=NA),
	namesData=list(CNOlist=NA,model=NA),
	resE=NULL){

 if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    } 





    if (is.null(resE)==TRUE){
        resE<-residualError(CNOlist)
    }

	
#Create a report directory and copy css and logos in there
	dir.create(directory)	
	cpfile<-dir(system.file("templates",package="CellNOptR"),full.names=TRUE)
	resultsdir<-file.path(getwd(),directory)
	file.copy(from=cpfile,to=resultsdir,overwrite=TRUE)
	
#copy all the graphs into the directory

	for(i in 1:length(namesFiles)){
	
		if(!is.na(namesFiles[[i]])){
			cpfile<-file.path(getwd(),namesFiles[[i]])
			file.copy(from=cpfile,to=resultsdir,overwrite=TRUE)
			unlink(cpfile)
			}
			
		}
		
#Produce the html page
	htmlfile=paste(resultsdir,'/CellNOptReport.html',sep="")
	cat(
		'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"http://www.w3.org/TR/html4/loose.dtd"> \n <html> ', 
		file = htmlfile, append=TRUE)
	cat(
		'\n <link rel="stylesheet" type="text/css" href="CNOR.css" /> \n ', 
		file = htmlfile, append=TRUE)
	cat(
		'<head> <title> CellNOptR Analysis Report </title> </head> \n ', 
		file = htmlfile, append=TRUE)
	cat(
		'<body> \n <table class="border"> \n <tr class="border top"> \n ', 
		file = htmlfile, append=TRUE)
	cat(
		'<td class="border corner"> &nbsp&nbsp&nbsp&nbsp </td>\n <td class="border top"> ', 
		file = htmlfile, append=TRUE)
	cat(
		'<h1>Report from CellNOptR &nbsp &nbsp </h1> </td> \n ', 
		file = htmlfile, append=TRUE)
	cat(
		'<td class="border corner2"> &nbsp&nbsp&nbsp&nbsp </td>\n </tr>\n <tr class="border"> \n <td class="border"></td>\n ', 
		file = htmlfile, append=TRUE)
	cat(
		'<td class="main">\n <table>\n <tr>\n ', 
		file = htmlfile, append=TRUE)
	cat(
		'<img src="Rlogo.png" width="70" height="70" ALIGN=RIGHT>&nbsp ', 
		file = htmlfile, append=TRUE)
	cat(
		'<img src= "phinou.jpeg" width="110" height="110" ALIGN=RIGHT>&nbsp ', 
		file = htmlfile, append=TRUE)
	cat(
		'<img src="logo_saezrodriguez_mid.gif" width="110" height="80" ALIGN=RIGHT>\n </tr>\n <tr>', 
		file = htmlfile, append=TRUE)
	cat(
		paste(
			'\n  <td class="tabs"> <h3><a href="',
			namesFiles$dataPlot,
			'" title="data" target=blank >Data plot</a></h3> </td> ',
			sep=""), 
		append = TRUE, file = htmlfile)
	cat(
		paste(
			'\n  <td class="tabs"> <h3><a href="',
			namesFiles$simResultsT1,
			'" title="Simulation results t1" target=blank >Simulation Results t1</a></h3> </td> ',
			sep=""), 
		append = TRUE, file = htmlfile)
	cat(
		paste(
			'\n  <td class="tabs"> <h3><a href="',
			namesFiles$evolFitT1,
			'" title="Evolution of fit t1" target=blank >Evolution of fit t1</a></h3> </td> ',
			sep=""), 
		append = TRUE, file = htmlfile)
		
	if(!is.na(namesFiles$evolFitT2)){
		cat(
			paste(
				'\n  <td class="tabs"> <h3><a href="',
				namesFiles$simResultsT2,
				'" title="Simulation results t2" target=blank >Simulation Results t2</a></h3> </td> ',
				sep=""), 
			append = TRUE, file = htmlfile)
		cat(
			paste(
				'\n  <td class="tabs"> <h3><a href="',
				namesFiles$evolFitT2,
				'" title="Evolution of fit t2" target=blank >Evolution of fit t2</a></h3> </td> ',
				sep=""), 
			append = TRUE, file = htmlfile)
		}
		
	cat('\n </tr> ', append = TRUE, file = htmlfile)	
	cat(
		paste('<tr> \n <br> General information: \n <UL> \n <LI>data: ',namesData$CNOlist,sep=""), 
		append = TRUE, file = htmlfile)
	cat(
		paste('\n <LI>time point(s): ',ifelse(is.na(namesFiles$evolFitT2),1,2),sep=""), 
		append = TRUE, file = htmlfile)
	cat(
		paste('\n <LI>Residual error: t1: ',resE["t1"],sep=""), 
		append = TRUE, file = htmlfile)
	nDataP<-sum(!is.na(CNOlist@signals[[2]]))
	scaledResE<-resE["t1"]/nDataP
	cat(
		paste(' (scaled: ',round(scaledResE,digits=4),' )',sep=""), 
		append = TRUE, file = htmlfile)
		
	if(!is.na(resE["t2"])) {
		cat(paste('; t2: ',resE["t2"],sep=""), append = TRUE, file = htmlfile)
		nDataP<-sum(!is.na(CNOlist@signals[[3]]))
		scaledResE<-resE["t2"]/nDataP
		cat(
			paste(' (scaled: ',round(scaledResE,digits=4),' )',sep=""), 
			append = TRUE, file = htmlfile)
		cat(
			paste('; t1andt2: ',resE["t1andt2"],sep=""), 
			append = TRUE, file = htmlfile)
		nDataP<-sum(!is.na(CNOlist@signals[[3]]))+sum(!is.na(CNOlist@signals[[2]]))
		scaledResE<-resE["t1andt2"]/nDataP
		cat(
			paste(' (scaled: ',round(scaledResE,digits=4),' )',sep=""), 
			append = TRUE, file = htmlfile)
		}
		
	cat(
		paste('\n <LI>previous knowledge network: ',namesData$model,sep=""), 
		append = TRUE, file = htmlfile)
	cat(
		paste(
			'\n <LI>PKN: ',
			dim(modelOriginal$interMat)[1],' species and ',
			dim(modelOriginal$interMat)[2],' interactions',sep=""), 
		append = TRUE, file = htmlfile)
	cat(
		paste(
			'\n <LI>Scaffold (compressed and expanded): ',
			dim(modelOpt$interMat)[1],' species and ',
			dim(modelOpt$interMat)[2],' interactions', '\n </UL>\n'), 
		append = TRUE, file = htmlfile)
	cat(
		paste(
			'<br> Optimisation t1: \n<UL> \n <LI>generations: ',
			dim(optimResT1$results)[1],' (best model obtained after ',
			optimResT1$results[dim(optimResT1$results)[1],"Stall_Generation"],') \n', sep=""),
		append = TRUE, file = htmlfile)
 	cat(
 		paste(
 			'<LI>best string: ',sum(optimResT1$bString),
 			' interactions, objective function = ',
 			optimResT1$results[dim(optimResT1$results)[1],
 			"Best_score_Gen"],'\n', sep=""),
 		append = TRUE, file = htmlfile)
 	cat(
 		paste(
 			'<LI>number of strings within the tolerance limits: ',
 			dim(optimResT1$stringsTol)[1],' \n</UL>\n', sep=""),
 		append = TRUE, file = htmlfile)
 		
	if(!is.na(namesFiles$evolFitT2)){
		cat(paste(
			'<br> Optimisation t2: \n<UL> \n <LI>generations: ',
			dim(optimResT2$results)[1],' (best model obtained after ',
			optimResT2$results[dim(optimResT2$results)[1],
			"Stall_Generation"],') \n', sep=""),append = TRUE, file = htmlfile)
 		cat(paste(
 			'<LI>best string: ',sum(optimResT2$bString),
 			' additional interactions, objective function = ',
 			optimResT2$results[dim(optimResT2$results)[1],"Best_score_Gen"],
 			'\n', sep=""),append = TRUE, file = htmlfile)
 		cat(paste(
 			'<LI>number of strings within the tolerance limits: ',
 			dim(optimResT2$stringsTol)[1],'\n</UL>\n', sep="")
 			,append = TRUE, file = htmlfile)
		}
		
 	cat(paste(
 		'<br> Scaffold network:\n <UL>\n <LI>cytoscape sif format: ',
 		namesFiles$scaffold,'\n',sep=""),append = TRUE, file = htmlfile)
 	cat(paste(
 		'<LI>graphviz dot format: ',namesFiles$scaffoldDot,'\n',sep=""),
 		append = TRUE, file = htmlfile)
 	cat(paste(
 		'<LI>edge attribute files: ',namesFiles$wscaffold,' and ',
 		namesFiles$tscaffold,'\n',sep=""),append = TRUE, file = htmlfile)
 	cat(paste(
 		'<LI>edge attributes respectively reflect the weight of the edge calculated as the frequency of the edge in the best solutions in the last generation of optimisation(s), and the time stamp on each edge (i.e. absent/present at t1/present at t2) \n</UL>\n',
 		sep=""),append = TRUE, file = htmlfile)
 	cat(paste(
 		'<br> Previous knowledge network:\n<UL> \n <LI>cytoscape sif format: ',
 		namesFiles$PKN,'\n',sep=""),append = TRUE, file = htmlfile)
 	cat(paste(
 		'<LI>graphviz dot format: ',namesFiles$PKNdot,'\n',sep=""),
 		append = TRUE, file = htmlfile)
 	cat(paste(
 		'<LI>edge attribute file: ',namesFiles$wPKN,' \n',sep=""),
 		append = TRUE, file = htmlfile)
 	cat(paste(
 		'<LI>node attribute file: ',namesFiles$nPKN,'\n',sep=""),
 		append = TRUE, file = htmlfile)
 	cat(paste(
 		'<LI>the edge attribute is a time stamp on each edge of the scaffold network (i.e. absent/present at t1/present at t2) mapped back to the original network\n',
 		sep=""),append = TRUE, file = htmlfile)
 	cat(paste(
 		'<LI>the node attribute reflects the information about which nodes are signals/inhibited/stimulated/non-controlable/non-observable\n </UL> \n',
 		sep=""),append = TRUE, file = htmlfile)
 	cat(paste(
 		'<br> <datestamp> ',date(),'</datestamp> ',sep=""),
 		append = TRUE, file = htmlfile)
 	cat(
 		'<hr/> \n </tr> \n </table> \n </td> \n <td class="border right"> </td> \n </tr> \n </table> \n </body> \n </html>',
 		append = TRUE, file = htmlfile)	 				
	}


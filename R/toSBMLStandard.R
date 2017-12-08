# Export a Boolean network <network> to an sbml-qual file <fileName>.
# This file can then be read in using other software that supports SBMLqual standards.
# The function also takes a bit string as input. 
#It cuts the model according to the values in bitstrings and write the new model object to SBMLqual.

library("stringi")
library("stringr")

toSBMLStandard <- function(network, file, bitString = c(rep(1,length(network$reacID))))
{
    
    network = cutModel(network, bitString)
    
    
    # generate a network identifier from the file name
    id <- sub(".sbml", "", basename(file), fixed=TRUE)
    id <- gsub("[^a-zA-Z0-9_]+","_",id)
    
    # open a string connection
    output <- NULL
    f <- textConnection("output",  encoding="UTF-8", open="w", local=TRUE)
    
    # write document header
    cat(file=f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    cat(file=f, "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version1/core\" level=\"3\" version=\"1\" xmlns:qual=\"http://www.sbml.org/sbml/level3/version1/qual/version1\" qual:required=\"true\">\n")
    cat(file=f, "\t<model id=\"", id , "\">\n", sep="")
    
    # write default compartment
    cat(file=f, "\t\t<listOfCompartments>\n")
    cat(file=f, "\t\t\t<compartment id=\"default\" constant=\"true\"/>\n")
    cat(file=f, "\t\t</listOfCompartments>\n")
    
    # write genes
    geneList <- NULL
    cat(file=f, "\t\t<qual:listOfQualitativeSpecies>\n")
    for (gene in network$namesSpecies)
    {
        Input_interaction <- which(network$interMat[gene,] == 1)
        Output_Interaction <- which(network$interMat[gene,] == -1)
        #if gene is has no input interaction and no output interaction, do not write it out to SBMLqual 
        if (length(Input_interaction) == 0 && length(Output_Interaction) == 0)next
        cat(file=f, "\t\t\t<qual:qualitativeSpecies qual:id =\"", gene, "\" qual:compartment=\"main\"",
            " qual:constant=\"false\"/>\n", sep = "")
        geneList <- c(geneList, gene)
        
    }
    cat(file=f, "\t\t</qual:listOfQualitativeSpecies>\n")
    
    # write transition functions
    cat(file=f, "\t\t<qual:listOfTransitions>\n")
    
    t_count = 1 
    
    for (gene in geneList)
    {
        inputs = c()
        signs = c()
        interactions <- which(network$interMat[gene,] == 1)
        t_name <- ""
        t_name = paste0("t", t_count, sep='')
        
        cat(file=f, "\t\t\t<qual:transition qual:id=\"", t_name , "\">\n",sep = "")
        cat(file=f,"\t\t\t\t <qual:listOfInputs> \n")
        
        for (i in interactions)
        {
            
            intName <- colnames(network$interMat)[i]
            if (grepl("+", intName,fixed = T) == TRUE) 
            {
                tmp = unlist(strsplit(intName,split = "="))[1] #take the inputs on the left hand side, removing the output gene
                LHS = unlist(strsplit(tmp,split = "+",fixed = T)) #list of inputs
                for (input in LHS){
                    sign <- "positive"
                    if(substr(input,1,1) == "!"){
                        input = stri_sub(input,2)
                        sign <- "negative"
                        
                    }
                    full_name <- ""
                    full_name <- paste0("theta_", t_name, "_", input, sep='')
                    signs <- c(signs, sign)
                    cat(file=f,"\t\t\t\t\t <qual:input  qual:id=\"", full_name , 
                        "\" qual:qualitativeSpecies=\"",input,"\" qual:transitionEffect=\"none\" qual:sign=\"", sign, "\" qual:thresholdLevel=\"1\"/>\n",sep = "")
                }
                andGate = paste(unlist(LHS), collapse='&')
                andGate = str_replace_all(andGate,'!','')
                if (length(inputs)!=0)
                    andGate = paste0('|', andGate)
                inputs = c(inputs, andGate)
                
            } 
            else 
            {
                sign <- "positive"
                LHS <- unlist(strsplit(intName,split = "="))[1] #single OR input
                if (substr(LHS,1,1) == "!"){
                    sign <- "negative"
                    LHS <- stri_sub(LHS,2)
                }
                
                full_name <- ""
                full_name <- paste0("theta_", t_name, "_",LHS,sep='')
                signs <- c(signs, sign)
                cat(file=f,"\t\t\t\t\t <qual:input  qual:id=\"", full_name , 
                    "\" qual:qualitativeSpecies=\"", LHS, "\" qual:transitionEffect=\"none\" qual:sign=\"", sign, "\" qual:thresholdLevel=\"1\"/>\n",sep = "")
                
                if (length(inputs)!=0)
                    orGate = paste0('|',LHS, collapse = '')
                else orGate = LHS
                inputs = c(inputs, orGate)
            }
        
            
        }
        cat(file=f,"\t\t\t\t </qual:listOfInputs>\n")
        cat(file=f,"\t\t\t\t <qual:listOfOutputs>\n")
        cat(file=f,"\t\t\t\t\t <qual:output qual:qualitativeSpecies=\"", gene ,"\" qual:transitionEffect=\"assignmentLevel\"/>\n",sep = "")
        cat(file=f,"\t\t\t\t </qual:listOfOutputs>\n")
        cat(file=f, "\t\t\t\t<qual:listOfFunctionTerms>\n")
        cat(file=f, "\t\t\t\t\t<qual:defaultTerm qual:resultLevel=\"0\"/>\n")
        cat(file=f, "\t\t\t\t\t<qual:functionTerm qual:resultLevel=\"1\">\n")
        cat(file=f, "\t\t\t\t\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n")
        cat(file=f, "\t\t\t\t\t\t\t<apply>\n")
        transition = paste(unlist(inputs), collapse='')
       
        #only OR transitions
        if (grepl("&", transition,fixed = T) == FALSE && transition!='')
        {
            transition = unlist(strsplit(transition,split = '|', fixed = TRUE))
            if (length(transition) > 1)
                cat(file=f, "\t\t\t\t\t\t\t\t<", 'or', "/>\n",sep = "")
            for (i in 1:length(transition)){
                cat(file=f, "\t\t\t\t\t\t\t\t<apply>\n")
                cat(file=f, "\t\t\t\t\t\t\t\t\t<eq/>\n")
                cat(file=f, "\t\t\t\t\t\t\t\t\t<ci>" , transition[i] , "</ci>\n")
                if (signs[i] == 'positive')
                    integer = 1
                else integer = 0
                cat(file=f, "\t\t\t\t\t\t\t\t\t<cn type='integer'>" , integer , "</cn>\n")
                cat(file=f, "\t\t\t\t\t\t\t\t</apply>\n")
            }
            
                    
                
        }
        #only AND transitions
        else if (grepl("&", transition,fixed = T) == TRUE && transition!='' && grepl("|", transition,fixed = T) == FALSE){
            transition = unlist(strsplit(transition,split = '&', fixed = TRUE))
            cat(file=f, "\t\t\t\t\t\t\t\t<", 'and', "/>\n",sep = "")
            for (i in 1:length(transition)){
                cat(file=f, "\t\t\t\t\t\t\t\t<apply>\n")
                cat(file=f, "\t\t\t\t\t\t\t\t\t<eq/>\n")
                cat(file=f, "\t\t\t\t\t\t\t\t\t<ci>" , transition[i] , "</ci>\n")
                if (signs[i] == 'positive')
                    integer = 1
                else integer = 0
                cat(file=f, "\t\t\t\t\t\t\t\t\t<cn type='integer'>" , integer , "</cn>\n")
                cat(file=f, "\t\t\t\t\t\t\t\t</apply>\n")
            
            }
        
        }
        else if (transition==''){
            #if gene has no input transition, do nothing 
        }
        #both And and OR gates 
        else 
        {
            count_sign = 1
            transition = unlist(strsplit(transition,split = '|', fixed = TRUE))
            cat(file=f, "\t\t\t\t\t\t\t\t<", 'or', "/>\n",sep = "")
            for (i in 1:length(transition)){
                cat(file=f, "\t\t\t\t\t\t\t\t<apply>\n")
                if (grepl("&", transition[i], fixed = T) == TRUE){
                    cat(file=f, "\t\t\t\t\t\t\t\t<", 'and', "/>\n",sep = "")
                    andTransition = unlist(strsplit(transition[i],split = '&', fixed = TRUE))
                    for (k in andTransition){
                        cat(file=f, "\t\t\t\t\t\t\t\t<apply>\n")
                        cat(file=f, "\t\t\t\t\t\t\t\t\t<eq/>\n")
                        cat(file=f, "\t\t\t\t\t\t\t\t\t<ci>" , k , "</ci>\n")
                        if (signs[count_sign] == 'positive')
                            integer = 1
                        else integer = 0
                        cat(file=f, "\t\t\t\t\t\t\t\t\t<cn type='integer'>" , integer , "</cn>\n")
                        cat(file=f, "\t\t\t\t\t\t\t\t</apply>\n")
                        count_sign = count_sign + 1 
                    }
                }
                else {
                    cat(file=f, "\t\t\t\t\t\t\t\t\t<eq/>\n")
                    cat(file=f, "\t\t\t\t\t\t\t\t\t<ci>" , transition[i] , "</ci>\n")
                    if (signs[i] == 'positive')
                        integer = 1
                    else integer = 0
                    cat(file=f, "\t\t\t\t\t\t\t\t\t<cn type='integer'>" , integer , "</cn>\n")
                    count_sign = count_sign + 1 
                }
                    
                cat(file=f, "\t\t\t\t\t\t\t\t</apply>\n")
            }
            
        }
        cat(file=f, "\t\t\t\t\t\t\t</apply>\n")
        cat(file=f, "\t\t\t\t\t\t</math>\n")
        cat(file=f, "\t\t\t\t\t</qual:functionTerm>\n")
        cat(file=f, "\t\t\t\t</qual:listOfFunctionTerms>\n")
        cat(file=f, "\t\t\t</qual:transition>\n")
        t_count = t_count + 1
    }
    
    # finish document
    cat(file=f, "\t\t</qual:listOfTransitions>\n")
    cat(file=f, "\t</model>\n")
    cat(file=f, "</sbml>\n")
    close(f)
    
    # open file and write the complete XML string
    f <- file(file, encoding="UTF-8", open="w")
    cat(file=f,output,sep="\n")
    close(f)  
    
}











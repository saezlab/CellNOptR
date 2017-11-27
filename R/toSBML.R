# Export a Boolean network <network> to an sbml-qual file <fileName>.

library("stringi")

writeLogic <- function(gene, inputs, t_name, t_count, logic, f){
    cat(file=f,"\t\t\t\t </qual:listOfInputs>\n")
    cat(file=f,"\t\t\t\t <qual:listOfOutputs>\n")
    cat(file=f,"\t\t\t\t\t <qual:output qual:qualitativeSpecies=\"", gene ,"\" qual:transitionEffect=\"assignmentLevel\"/>\n",sep = "")
    cat(file=f,"\t\t\t\t </qual:listOfOutputs>\n")
    cat(file=f, "\t\t\t\t<qual:listOfFunctionTerms>\n")
    cat(file=f, "\t\t\t\t\t<qual:defaultTerm qual:resultLevel=\"0\"/>\n")
    cat(file=f, "\t\t\t\t\t<qual:functionTerm qual:resultLevel=\"1\">\n")
    cat(file=f, "\t\t\t\t\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n")
    cat(file=f, "\t\t\t\t\t\t\t<apply>\n")
    
    if (length(inputs)!=0){
        if (length(inputs)>1)
            cat(file=f, "\t\t\t\t\t\t\t\t<",logic,"/>\n",sep = "")
        for (i in inputs){
            full_name <- ""
            full_name <- paste0("theta_", t_name, "_", i ,sep='')
            cat(file=f, "\t\t\t\t\t\t\t\t<apply>\n")
            cat(file=f, "\t\t\t\t\t\t\t\t\t<geq/>\n")
            cat(file=f, "\t\t\t\t\t\t\t\t\t<ci>" , i , "</ci>\n")
            cat(file=f, "\t\t\t\t\t\t\t\t\t<ci>" , full_name , "</ci>\n")
            cat(file=f, "\t\t\t\t\t\t\t\t</apply>\n")
        }
    }
    cat(file=f, "\t\t\t\t\t\t\t</apply>\n")
    cat(file=f, "\t\t\t\t\t\t</math>\n")
    cat(file=f, "\t\t\t\t\t</qual:functionTerm>\n")
    cat(file=f, "\t\t\t\t</qual:listOfFunctionTerms>\n")
    cat(file=f, "\t\t\t</qual:transition>\n")
    if (logic == "and") 
        t_count = t_count + 1
    return(t_count)
    
    
}

toSBML <- function(network, file)
{
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
    cat(file=f, "\t\t<qual:listOfQualitativeSpecies>\n")
    for (gene in network$namesSpecies)
    {

        cat(file=f, "\t\t\t<qual:qualitativeSpecies qual:id =\"", gene, "\" qual:compartment=\"main\"",
            " qual:constant=\"false\"/>\n", sep = "")
    
    }
    cat(file=f, "\t\t</qual:listOfQualitativeSpecies>\n")
    
    # write transition functions
    cat(file=f, "\t\t<qual:listOfTransitions>\n")
    
    t_count = 1 
    
    for (gene in network$namesSpecies)
    {
        #keep track of all interaction for a gene
        all_interactions <- which(network$interMat[gene,]==1)
        or_interaction = c()
        and_interaction = c()
        #write all OR interactions 
        t_name <- ""
        t_name = paste0("t", t_count, sep='')
        cat(file=f, "\t\t\t<qual:transition qual:id=\"", t_name , "\">\n",sep = "")
        cat(file=f,"\t\t\t\t <qual:listOfInputs> \n")
        l = list()
        
        for (i in all_interactions){
            int_name <- colnames(network$interMat)[i]
            if (grepl("+", int_name,fixed = T) == TRUE) {
                and_interaction = c(and_interaction,int_name)
                next
            }
            else{
                
                sign <- "positive"
                LHS <- unlist(strsplit(int_name,split = "="))[1] #input
                if (substr(LHS,1,1) == "!"){
                    sign <- "negative"
                    LHS <- stri_sub(LHS,2)
                }
                or_interaction = c(or_interaction, LHS)
                RHS <- unlist(strsplit(int_name,split = "="))[2] #gene -- output
                full_name <- ""
                full_name <- paste0("theta_", t_name, "_",LHS,sep='')
                cat(file=f,"\t\t\t\t\t <qual:input  qual:id=\"", full_name , 
                    "\" qual:qualitativeSpecies=\"", LHS, "\" qual:transitionEffect=\"none\" qual:sign=\"", sign, "\" qual:thresholdLevel=\"1\"/>\n",sep = "")
            }
        }
        #write the logical headers for the OR transitions
        
        writeLogic(gene,or_interaction,t_name, t_count, "or",f)
       
        #if there are AND interactions, create new transition
        if (length(and_interaction)==0)
            t_count = t_count + 1
        for (i in and_interaction){
            t_count = t_count + 1
            and_list = c()
            t_name <- ""
            t_name = paste0("t", t_count, sep='')
            cat(file=f, "\t\t\t<qual:transition qual:id=\"", t_name ,  "\">\n",sep = "")
            cat(file=f,"\t\t\t\t <qual:listOfInputs> \n")
            tmp = unlist(strsplit(i,split = "="))[1] #input
            LHS = unlist(strsplit(tmp,split = "+",fixed = T))
            for (input in LHS){
                sign <- "positive"
                if(substr(input,1,1) == "!"){
                    input = stri_sub(input,2)
                    sign <- "negative"
                }
                and_list = c(and_list,input)
                full_name <- ""
                full_name <- paste0("theta_", t_name, "_", input, sep='')
                cat(file=f,"\t\t\t\t\t <qual:input  qual:id=\"", full_name , 
                    "\" qual:qualitativeSpecies=\"",input,"\" qual:transitionEffect=\"none\" qual:sign=\"", sign, "\" qual:thresholdLevel=\"1\"/>\n",sep = "")
            }
            t_count = writeLogic(gene,and_list,t_name, t_count, "and",f)
            
        }
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


    








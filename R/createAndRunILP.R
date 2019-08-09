#
#  This file is part of the CNO software
#
#  Copyright (c) 2019 - SaezLab - Heidelberg University
#
#  File author(s): CNO developers
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################

createAndRunILP <- function(model = model, 
                            midas = midas, 
                            cnolist = cnolist, 
                            accountForModelSize = accountForModelSize, 
                            sizeFac = sizeFac, 
                            source_path = source_path, 
                            mipGap = mipGap, 
                            relGap = relGap, 
                            timelimit = timelimit,
                            cplexPath = cplexPath, 
                            method = method, 
                            numSolutions = numSolutions, 
                            limitPop = limitPop, 
                            poolIntensity = poolIntensity, 
                            poolReplace = poolReplace){
  
  # Creating auxilliary variables
  time_start_writing = Sys.time()
  md <- midas
  TimeOfEvaluation = max(md$dataMatrix[,md$DAcol[1]])
  
  # Creating the vector of edge identifiers
  create_y_vector_from_preprocessed_model <- function(model){
    
    y_vector <- model$reacID
    return(y_vector)
    
  }
  
  # Creating set of variables to inhibitory effects of each interactions
  create_inhibitor_set_for_y_i <- function(model, y_i){
    
    inhibitor_set = c()
    inhibitor_set = names(which(model$notMat[,which(model$reacID == y_i)] =="1"))
    return(inhibitor_set)
    
  }
  
  # Creating set of variables to the reactants of each interactions
  create_reactant_set_for_y_i <- function (model, y_i, inhibitor_set){
    
    reactant_set = c()
    reactant_set = names(which(model$interMat[,which(model$reacID == y_i)] == "-1"))
    remove <- inhibitor_set
    reactant_set = reactant_set[!reactant_set %in% remove] #removes inhibitors from the reactant set (are showing up in the interMat, so they are originally considered to be reactants)
    return(reactant_set)
    
  }
  
  # Creating set of variables to products of each interactions
  create_product_set_for_y_i <- function(model, y_i){
    
    product_set = c()
    product_set = names(which(model$interMat[,which(model$reacID == y_i)] == "1"))
    return(product_set)
    
  }
  
  # Combining all interaction set of variables
  create_reaction_for_y_i <- function(model, y_i){
    
    inhibitor_set = create_inhibitor_set_for_y_i(model, y_i) 
    reactant_set = create_reactant_set_for_y_i(model, y_i, inhibitor_set) # also removes inhibitors that are present in the inhibitor set of the respective reaction, that's why the second argument is necessary
    product_set = create_product_set_for_y_i(model, y_i)
    reaction_set = list(reactant_set, inhibitor_set, product_set)
    return(reaction_set)
    
  }
  
  # Creating all the sets
  create_all_sets <- function(model){ # sets: [[1]]: reactant set; [[2]]: inhibitor set; [[3]]: product set
    
    y_vector = create_y_vector_from_preprocessed_model(model)
    all_reaction_sets = data.frame(reac_name = c("reactant_set", "inhibitor_set", "product_set"))
    for(i in 1:length(y_vector)){
      y_i = y_vector[i]
      reaction_for_y_i = create_reaction_for_y_i(model, y_i)
      all_reaction_sets$y_1 <- reaction_for_y_i
      colnames(all_reaction_sets)[i+1] <- y_vector[i]
    }
    all_reaction_sets[[1]] = NULL
    return(all_reaction_sets)
    
  }
  
  # Extracting only the measurements from the MIDAS file
  extractMeasurementPartFromMidas <- function(midas, TimeOfEvaluation){
    
    md = midas$dataMatrix
    # md_out = md[which(md[midas$DAcol[1]] ==TimeOfEvaluation),midas$DVcol] # Accept also single measurement - added 13.02.18
    if (length(midas$DAcol)>1) {
      md_out <- md[which(md[midas$DAcol[1]] ==TimeOfEvaluation),midas$DVcol]
    } else { 
      md_out<-matrix(NA,length(which(md[midas$DAcol[1]] ==TimeOfEvaluation)),1);
      md_out[,1]<-md[which(md[midas$DAcol[1]] ==TimeOfEvaluation),midas$DVcol];
      colnames(md_out) <- colnames(md)[midas$DAcol]
    }
    return(md_out)
    
  }
  
  # Extracting only the treatment part from the MIDAS file
  extractTreatmentPartFromMidas <- function(midas, TimeOfEvaluation){
    
    md = midas$dataMatrix
    md_out = md[which(md[midas$DAcol[1]] ==TimeOfEvaluation),midas$TRcol]
    return(md_out)
    
  }
  
  # Means of measurements
  meansOfMeasurements <- function(midasExperimentPart_at_t0, midasExperimentPart, midasTreatmentPart){
    
    ValueMatrix = midasExperimentPart_at_t0
    if (!is.null(dim(ValueMatrix))) { # In case measurement at t0 is not a matrix (Added 13.02.18)
      if(all(midasTreatmentPart[1,]==0)){ #That means that the first experiment of the Treatment section has no stimuli --> measurements for t_1 for the first have to be taken into the mean (this is the "+1" in "N+1")
        ValueMatrix = rbind(ValueMatrix, midasExperimentPart[1,])
      }
      out = colMeans(ValueMatrix, na.rm = TRUE)
    } else {
      out = ValueMatrix
    }
    return(out)
    
  }
  
  # Creating and running the ILP problem
  reaction_sets = create_all_sets(model)
  midasExperimentPart = extractMeasurementPartFromMidas(md, TimeOfEvaluation)
  midasExperimentPart_at_t0 = extractMeasurementPartFromMidas(md,0)
  midasTreatmentPart = extractTreatmentPartFromMidas(md, TimeOfEvaluation)
  meansOfMeasurements_at_t0 = meansOfMeasurements(midasExperimentPart_at_t0, midasExperimentPart, midasTreatmentPart)
  numberOfExperiments = if(!is.null(dim(midasExperimentPart)[1])){dim(midasExperimentPart)[1]}else{1}
  y_vector = create_y_vector_from_preprocessed_model(model)
  binary_variables = create_binaries(model, md, numberOfExperiments, y_vector)
  objectiveFunction = writeObjectiveFunction(model, midasExperimentPart, y_vector = y_vector, accountForModelSize, sizeFac, meansOfMeasurements_at_t0, method = method)
  constraints = write_constraints(model, midasExperimentPart, midasTreatmentPart, reaction_sets, y_vector, md)
  bounds = write_bounds(model, midasTreatmentPart,y_vector, binary_variables)
  writeFile(objectiveFunction, constraints, bounds, binary_variables, cplexPath = cplexPath)
  time_cplex_start = Sys.time()
  invokeCPLEX(inputFileName = "testFile.lp", outputFileName = "cplexSolution.txt", mipGap=mipGap, relGap = relGap, timelimit=timelimit, cplexPath=cplexPath,
              numSolutions=numSolutions, limitPop = limitPop, poolIntensity = poolIntensity, poolReplace = poolReplace)
  time_cplex_end = Sys.time()
  bitstringILPAll = createILPBitstringAll("cplexSolution.txt", y_vector, binary_variables) #if interested in the solution for all binaries, use "readOutBinaries()" function
  
  time_end_writing = Sys.time()
  total_time = difftime(time_end_writing, time_start_writing, units = "secs")
  
  time_cplex_only = difftime(time_cplex_end, time_cplex_start, units = "secs")
  bScore = c()
  for(ii in 1:length(bitstringILPAll)){
    bScore = computeScoreT1(CNOlist = cnolist, model = model, bString = bitstringILPAll[[ii]])
  }
  
  return(list(bitstringILP = bitstringILPAll, bScore = bScore, time_cplex_only = time_cplex_only, total_time = total_time))
  
}

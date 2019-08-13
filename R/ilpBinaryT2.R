#
#  This file is part of the CNO software
#
#  Copyright (c) 2019 - SaezLab - Heidelberg Universit
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

ilpBinaryT2 <- function(cnolist = cnolist,
                        model = model, 
                        sizeFac = 0.0001, 
                        mipGap = 0, 
                        relGap = 0, 
                        timelimit = 3600, 
                        cplexPath = "~/Documents/cplex", 
                        method = "quadratic",
                        numSolutions = 100, 
                        limitPop = 500, 
                        poolIntensity = 0, 
                        poolReplace = 2){
  
  # Defining funtion that writes the objective function of the ILP problem
  writeObjectiveFunction <- function(model, 
                                     midasExperimentPart, 
                                     y_vector=y_vector,
                                     accountForModelSize = TRUE, 
                                     sizeFac = .000001, 
                                     meansOfMeasurements_at_t0, 
                                     method = "quadratic" ){ # function returns 2 values in a list: 1.) part of the objetive that only consists of absolute values, 2.) constant offset.
    
    # creating auxilliary objects
    reducedMidas = midasExperimentPart
    whichObjective = method
    objective1 = c()
    
    # need link between column header and species in the binary_variables vector
    link_midas_colnames_to_namesSpecies <- function(model, reducedMidas){
      
      linker = c()
      if (!is.null(colnames(reducedMidas))) { # if there is no column name (single experiment) - added 13.02.18
        for(i in 1:length(colnames(reducedMidas))){
          if(str_sub(colnames(reducedMidas)[i],-1)=="i"){# for the inhibitors in the treatement section of the  midas file
            linker[i] = which(model$namesSpecies==substr(strsplit(colnames(reducedMidas)[i],":")[[1]][2],1,nchar(strsplit(colnames(reducedMidas)[i],":")[[1]][2])-1))
          } else{  #for everything else
            linker[i] = which(model$namesSpecies==strsplit(colnames(reducedMidas)[i],":")[[1]][2])
          }
        }
      } else {
        linker <- 1:length(model$namesSpecies)
      }
      return(linker)
      
    }
    
    # writing the objective function
    linker_vector = link_midas_colnames_to_namesSpecies(model, reducedMidas)
    weight_of_unstimulated_experiments = if(!is.null(dim(reducedMidas)[1])){dim(reducedMidas)[1]+1}else{1} # this is N+1 / if clause added 13.02.18 
    constant = 0
    factors = "obj:\t"
    if(whichObjective == "linear"){
      for(i in 1:dim(reducedMidas)[1]){
        
        for(j in 1:dim(reducedMidas)[2]){
          if(!is.na(reducedMidas[i,j])){ #if measurements are NA (= not measured) they cannot appear in the objective function
            #term1:
            if(i ==1){
              term1 = weight_of_unstimulated_experiments * meansOfMeasurements_at_t0[j]
            }else{
              term1 = reducedMidas[i,j]
            }
            #term2:
            if(i ==1){
              term2 = weight_of_unstimulated_experiments * (1-2*meansOfMeasurements_at_t0[j])
            }else{
              term2 = 1-2*reducedMidas[i,j]
            }
            #term3:
            #here, the number of the binary is needed - procedure: use linker_vector to check, which number the respective species has in the namesSpecies Vector, then add 
            #(1.) the number of reactions (=offset), 
            #(2.) the number of experiments that have been evaluated before 
            #(3.) the index of the species evaluated in this experiment, 
            #(4.) (-1) since the binaries start with index 0 instead of 1.
            #dim(reducedMidasFile)[1]
            binary_index = length(model$reacID) + (i-1)*length(model$namesSpecies) + linker_vector[j]  # -1 is used since the binaries in Mitsos code start with xb0 - can be taken out after prove of principle
            term3 = paste0("xb", binary_index)#, "    ", binary_variables[[2]][binary_index])
            constant = constant + term1
            
            #print(paste0("Variable: xb", binary_index, " term1: ", term1, " term2: ", term2, " term3: ", term3, " meansOfMeasurements_at_t0: ", meansOfMeasurements_at_t0[j]))
            
            if(term2<0){
              factors  = paste(factors,term2,term3)
            }else{
              factors = paste(factors, "+", term2, term3)
            }
            
          }
          
          
        }
        if(accountForModelSize){
          for(m in 1:length(y_vector)){
            factors = paste0(factors, " + ", sizeFac, " xb", m )
          }
        }
        
      }
      # here, the model size is accounted for:
      for(m in 1:length(y_vector)){
        factors = paste0(factors, " + ", sizeFac, " xb", m )
      }
      
      out = list(constant, factors) #see comment in the functon definition above
      
      
    } else if(whichObjective =="quadratic")
    {
      for(i in 1:dim(reducedMidas)[1]){
        # for boolean predictions, the square difference between measurement and prediction can be calculated as
        # (x-x_m)^2 = x -   2*x*x_m +   x_m*x_m
        #           term_1    term2       term3
        for(j in 1:dim(reducedMidas)[2]){
          
          if(!is.na(reducedMidas[i, j])){
            
            binary_index = length(model$reacID) + (i-1)*length(model$namesSpecies) + linker_vector[j]  # -1 is used since the binaries in Mitsos code start with xb0 - can be taken out after prove of principle
            #real term1:
            term1 = paste0("xb", binary_index)
            #term2:
            if(reducedMidas[i, j]>=0){
              term2 = paste0(2*abs(reducedMidas[i,j]), " xb", binary_index)
            } else {
              term2 = paste0(0, " xb", binary_index)
            }
            #real term3:
            term3 = reducedMidas[i,j] *reducedMidas[i,j]
            
            factors  = paste(factors,"+", term1,"-", term2)#,"+", term3)
            
          }
          
          
        }
        
        
        out = list(constant, factors) #see comment in the functon definition above
      }
    }
    
    return(out)
    
  }
  
  # Function that creates the binary variables of the ILP problem
  create_binaries <- function(model, 
                              midas, 
                              numberOfExperiments, 
                              y_vector){
    
    # creating auxilliary variables
    numbers = c()
    variables = c()
    identifiers = c()
    
    # function for creating binary variables for the edges
    create_binary_variables_for_y_vector <- function(y_vector){
      
      variables = c()
      for(i in 1:length(y_vector)){
        variables[i] = paste0("y_",i)
      }
      binary_variables_list = list(1:(length(y_vector)),variables, paste("reaction",y_vector))
      return(binary_variables_list)
      
    }
    
    # function for creating binary variables for the species at each experimental condition separately
    create_binary_variables_for_all_species_in_experiment_k_i <- function(model, 
                                                                          numberOfExperiment){ # k means experiment index
      
      number_of_species = length(model$namesSpecies)
      identifier = c()
      variables = c()
      for(i in 1:number_of_species){
        identifier[i] = paste("species", model$namesSpecies[i], "in experiment k =", numberOfExperiment)
        variables[i] = paste0("x_",i,"^",numberOfExperiment)
      }
      binary_variables_list = list(0:length(model$namesSpecies),variables, identifier)
      return(binary_variables_list)
      
    }
    
    # function for creating binary variables for the interaction state at each experimental condition separately
    create_binary_variables_for_all_reactions_z_i_in_experiment_k_i <- function(model, 
                                                                                numberOfExperiment){ # k means experiment index
      number_of_reactions = length(model$reacID)
      identifier = c()
      variables = c()
      for(i in 1:number_of_reactions){
        identifier[i] = paste("reaction", model$reacID[i], "in experiment k =", numberOfExperiment)
        variables[i] = paste0("z_",i,"^",numberOfExperiment)
      }
      binary_variables_list = list(0:length(model$reacID),variables, identifier)
      return(binary_variables_list)
    }
    
    # Creating binary variables for the edges
    y_binaries = create_binary_variables_for_y_vector(y_vector)
    numbers = y_binaries[[1]]
    variables = y_binaries[[2]]
    identifiers = y_binaries[[3]]
    
    # creating binary variables for the species
    for(i in 1:numberOfExperiments){
      binList = create_binary_variables_for_all_species_in_experiment_k_i(model, i)
      numbers= append(numbers, binList[[1]])
      variables = append(variables, binList[[2]])
      identifiers = append(identifiers, binList[[3]])
    }
    for(i in 1:numberOfExperiments){
      binList = create_binary_variables_for_all_reactions_z_i_in_experiment_k_i(model, i)
      numbers= append(numbers, binList[[1]])
      variables = append(variables, binList[[2]])
      identifiers = append(identifiers, binList[[3]])
    }
    
    # combining and returning all the binary variables
    binaries = list(c(paste0("xb",1:length(variables))), variables, identifiers)
    return(binaries)
    
  }
  
  # Function that write the bounds of th ILP variables
  write_bounds <- function(model,
                           midasTreatmentPart,
                           y_vector,
                           binary_variables){
    
    # package stringr is needed in order to run this function. otherwise problems in the str_sub() command
    # this function returns the bounds vector which contains the boundaries for each binary varible.
    # procedure: create a vector that has the default binary bounds 0 <= xb_i <= 1
    # then, overwrite the elements of the vector that are fixed to a number (e.g. Treatment 1/0 or inhibitor 1) with the respective value.
    
    # writing the bounds for the edge variables (0 or 1)
    writeBoundsForEdges_y_i<- function(y_vector){
      y_bound_vector = c()
      for(i in 1:length(y_vector)){
        y_bound_vector[i] = paste0("0 <= xb", i, " <= 1")
      }
      return(y_bound_vector)
    }
    
    # linking species names in the midas with the ones in the model
    link_midas_colnames_to_namesSpecies <- function(model, reducedMidas){
      
      linker = c()
      if (!is.null(colnames(reducedMidas))) { # if there is no column name (single experiment) - added 13.02.18
        for(i in 1:length(colnames(reducedMidas))){
          if(str_sub(colnames(reducedMidas)[i],-1)=="i"){# for the inhibitors in the treatement section of the  midas file
            linker[i] = which(model$namesSpecies==substr(strsplit(colnames(reducedMidas)[i],":")[[1]][2],1,nchar(strsplit(colnames(reducedMidas)[i],":")[[1]][2])-1))
          } else{  #for everything else
            linker[i] = which(model$namesSpecies==strsplit(colnames(reducedMidas)[i],":")[[1]][2])
          }
        }
      } else {
        linker <- 1:length(model$namesSpecies)
      }
      return(linker)
      
    }
    
    # writing of the boundaries
    y_bounds_vector = writeBoundsForEdges_y_i(y_vector)
    linkerVector = link_midas_colnames_to_namesSpecies(model, midasTreatmentPart)
    bounds = c(paste0("0<= xb",1:(length(binary_variables[[1]])), " <= 1"))#, "   ", binary_variables[[2]][1:length(binary_variables[[3]])]))
    treatmentMatrix <- midasTreatmentPart
    for(i in 1:dim(midasTreatmentPart)[1]){
      for(j in 1:dim(midasTreatmentPart)[2]){
        if(midasTreatmentPart[i,j] == 1 && !(str_sub(colnames(midasTreatmentPart)[j], -1) == "i")){
          bounds[length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j]] = paste0("xb", length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j], " = 1")
        } else if(midasTreatmentPart[i,j] == 0 && !(str_sub(colnames(midasTreatmentPart)[j], -1) == "i")){
          bounds[length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j]] = paste0("xb", length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j], " = 0")
          
        }
        else if(midasTreatmentPart[i,j] == 1 && str_sub(colnames(midasTreatmentPart)[j], -1) == "i"){
          bounds[length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j]] = paste0("xb", length(y_vector) + (i-1)*length(model$namesSpecies)+linkerVector[j], " = 0")
        } else{
          #nothing. stays at it was.
          #bounds[length(y_vector) + (i-1)*length(model$namesSpecies) + linkerVector[j]-1] =
        }
      }
    }
    
    return(bounds)
  }
  
  # Function that writes the constraints of the ILP formulation
  write_constraints <- function(model,
                                midasExperimentPart,
                                midasTreatmentPart,
                                reaction_sets,
                                y_vector, 
                                midas, 
                                binary_variables){
    
    # constraint numbers refer to Mitsos et al., PloS Comp Biol. 2009
    # in th .lp file, the constraint numbers are written for debugging
    
    # constraint 3
    writeConstraint3 <- function(model, midasExperimentPart){#, numberOfExperiment, binary_variables){ # constraints are numbered according to the nomenclature in Mitsos 2009, PLoS Comp. Biol.
      
      reducedMidas = midasExperimentPart
      constraints = c()
      counter = 1
      for(i in 1:dim(reducedMidas)[1]){
        for(j in 1:length(model$reacID)){
          newConstraint = paste0("c", counter, ":\t", "xb", (length(model$reacID)+dim(reducedMidas)[1]*length(model$namesSpecies)+length(model$reacID)*(i-1)+j), " - xb" , j, " <= 0")#,"   ", binary_variables[[2]][(length(model$reacID)+dim(reducedMidas)[1]*length(model$namesSpecies)+length(model$reacID)*(i-1)+j)], " - ",  binary_variables[[2]][j] ,"  <=  0")
          constraints[counter] =newConstraint
          counter= counter+1
        }
      }
      
      return(constraints)
      
    }
    constraintsVector = writeConstraint3(model, midasExperimentPart)
    
    # constraints 4, 5 & 7
    reducedMidas = midasTreatmentPart
    counter = length(constraintsVector)+1
    for(i in 1:dim(reducedMidas)[1]){
      for(j in 1:length(model$reacID)){
        y_index = j
        if(length(reaction_sets[[j]][[1]])!=0){ #constraint 4
          for(k in 1:length(reaction_sets[[j]][[1]])){ #[[1]] for reactant set, [[2]] for inhibitor set, [[3]] for product set
            z_index = length(model$reacID) + length(model$namesSpecies)*dim(reducedMidas)[1] + (i-1)*length(model$reacID) + j 
            x_index = length(model$reacID) + (i-1)*length(model$namesSpecies) + which(model$namesSpecies==reaction_sets[[j]][[1]][k])
            newConstraint = paste0("c", counter, ":\t", "xb",z_index, " - xb",x_index, " <= 0 \t \t \\ c4")#, "    ", binary_variables[[2]][which(binary_variables[[1]] ==paste0("xb", z_index))], " - ", binary_variables[[2]][which(binary_variables[[1]] ==paste0("xb", x_index))], " <= 0")
            constraintsVector[counter] = newConstraint
            counter= counter+1
          }
        }
        if(length(reaction_sets[[j]][[2]])!=0){ #constraint 5
          for(k in 1:length(reaction_sets[[j]][[2]])){ #[[1]] for reactant set, [[2]] for inhibitor set, [[3]] for product set
            z_index = length(model$reacID) + length(model$namesSpecies)*dim(reducedMidas)[1] + (i-1)*length(model$reacID) + j 
            x_index = length(model$reacID) + (i-1)*length(model$namesSpecies) + which(model$namesSpecies==reaction_sets[[j]][[2]][k])
            newConstraint = paste0("c", counter, ":\t", "xb",z_index, " + xb",x_index, " <= 1\t \t \\ c5")#, "   ", binary_variables[[2]][which(binary_variables[[1]] ==paste0("xb", z_index))], " + ", binary_variables[[2]][which(binary_variables[[1]] ==paste0("xb", x_index))], " <= 1")
            constraintsVector[counter] = newConstraint
            counter= counter+1
          }
        }
        if(length(reaction_sets[[j]][[3]])!=0){ #constraint 7
          for(k in 1:length(reaction_sets[[j]][[3]])){ #[[1]] for reactant set, [[2]] for inhibitor set, [[3]] for product set
            speciesIsInhibitedInThisExperiment = FALSE
            if(length(which(colnames(reducedMidas) ==paste0("TR:",reaction_sets[[j]][[3]][k], "i")))==1){ # this workaround has to be done since the "which" command returns an "empty" which cannot be dealt with. 
              speciesIsInhibitedInThisExperiment = (reducedMidas[i,which(colnames(reducedMidas) ==paste0("TR:",reaction_sets[[j]][[3]][k], "i"))] == 1)
            }
            
            if(!speciesIsInhibitedInThisExperiment){ # if the species is perturbed to zero (e.g. inhibitors), this equation is obsolete und leads to a wrong constraint set. Cf. Written documentation (Hermann's Master thesis)
              z_index = length(model$reacID) + length(model$namesSpecies)*dim(reducedMidas)[1] + (i-1)*length(model$reacID) + j
              x_index = length(model$reacID) + (i-1)*length(model$namesSpecies) + which(model$namesSpecies==reaction_sets[[j]][[3]][k])
              newConstraint = paste0("c", counter, ":\t", "xb",z_index, " - xb",x_index, " <= 0 \t \t \\ c7")#, "   ", binary_variables[[2]][which(binary_variables[[1]] ==paste0("xb", z_index))], " - ", binary_variables[[2]][which(binary_variables[[1]] ==paste0("xb", x_index))], " <= 0")
              constraintsVector[counter] = newConstraint
              counter= counter+1
            }
            
          }
        }
        
        
        #constraint 6
        z_index = length(model$reacID) + length(model$namesSpecies)*dim(reducedMidas)[1]+ length(model$reacID)*(i-1) + j
        x_reactant_index = c()
        x_inhibitor_index = c()
        if(length(reaction_sets[[j]][[1]])!=0){# reactant set is not empty
          for(k in 1:length(reaction_sets[[j]][[1]])){
            x_reactant_index = append(x_reactant_index, length(model$reacID) + (i-1)*length(model$namesSpecies) + which(model$namesSpecies==reaction_sets[[j]][[1]][k]))
          }
        }
        
        if(length(reaction_sets[[j]][[2]])!=0){# inhibitors set is not empty
          for(k in 1:length(reaction_sets[[j]][[2]])){
            x_inhibitor_index  = append(x_inhibitor_index, length(model$reacID) + (i-1)*length(model$namesSpecies) + which(model$namesSpecies==reaction_sets[[j]][[2]][k]))
          }
        }
        
        newConstraint = paste0("xb", y_index, " - ", "xb", z_index)
        reactants_sum = paste0(" <= ", length(x_reactant_index)) # if there is no element in the x_inhibitor_index vector (i.e. length=0), this declaration persists and results in " <= 0". otherwise, it is overwritten by the for loop below. 
        inhibitors_sum = ""
        if(!is.null(x_reactant_index)){
          binary_sum = ""
          for(l in 1:length(x_reactant_index)){
            binary_sum = paste0(binary_sum, " + xb", x_reactant_index[l])
          }
          reactants_sum = paste0( binary_sum, " <= " , length(x_reactant_index))
        }
        if(!is.null(x_inhibitor_index)){
          binary_sum = ""
          for(l in 1:length(x_inhibitor_index)){
            binary_sum = paste0(binary_sum, " - xb", x_inhibitor_index[l])
          }
          inhibitors_sum = binary_sum
        }
        
        newConstraint = paste0("c", counter,":\t", newConstraint, inhibitors_sum, reactants_sum, " \t \t \\ c6")
        constraintsVector[counter] = newConstraint
        counter= counter+1
        # --- End Constraint 6
      }
      
      #new constraint 8
      for(j in 1:length(model$namesSpecies)){
        flag_vector = c()
        newConstraint = paste0("xb", length(y_vector)+(i-1)*length(model$namesSpecies)+j)
        for(k in 1:length(model$reacID)){
          
          if(reaction_sets[[k]][3] == model$namesSpecies[j]){
            flag_vector = c(flag_vector,k)
            #}
          }
        }
        added_reactions = ""
        if(!is.null(flag_vector)){
          for(k in 1:length(flag_vector)){
            added_reactions = paste0(added_reactions,"- xb", length(y_vector) + (length(model$namesSpecies))*dim(reducedMidas)[1] + (i-1)*length(model$reacID) + flag_vector[k])
          }
          newConstraint = paste0("c", counter, ": ", newConstraint, added_reactions, " <= 0" ," \t \t \\ c 8 new")
          constraintsVector[counter] = newConstraint
          counter= counter+1  
        }
      }#end new 8
    }
    
    return(constraintsVector)
  }
  
  # Function that writes the ILP problem
  writeFile <- function(objectiveFunction,
                        constraints,
                        bounds,
                        binaries,
                        cplexPath){
    
    data = paste0(gsub(pattern = "cplex", replacement = "", x = cplexPath), "/testFile.lp")
    write("enter Problem", data)
    write("", data, append = TRUE)
    write("Minimize", data, append = TRUE)
    
    write(objectiveFunction[[2]], data, append = TRUE)
    
    write("Subject To", data, append = TRUE)
    write(constraints, data, append = TRUE)
    
    write("Bounds", data, append = TRUE)
    write(bounds, data, append = TRUE)
    
    write("Binaries",data, append = TRUE)
    write(binaries[[1]], data, append = TRUE)
    
    write("End", data, append = TRUE)
    write("optimize Problem", data, append = TRUE)
    write("write solution_cplex_primaryObjective.txt", data, append = TRUE)
    write("sol", data, append = TRUE)
    write("y", data, append = TRUE)
    write("quit", data, append = TRUE)
    
  }
  
  # Function that executes the defined CPLEX problem
  invokeCPLEX <- function(inputFileName, 
                          outputFileName, 
                          mipGap=mipGap, 
                          relGap = relGap, 
                          timelimit=timelimit, 
                          cplexPath = cplexPath, 
                          numSolutions = numSolutions, 
                          limitPop = limitPop, 
                          poolIntensity = poolIntensity, 
                          poolReplace = poolReplace){
    
    # Function that writes the controll file for CPLEX. If needed, here additional options such as maximal computation time can be introduced
    writeCplexCommandFile <- function(commandFileName, 
                                      inputFileName, 
                                      outputFileName,
                                      mipGap=mipGap,
                                      relGap=relGap,
                                      timelimit=timelimit, 
                                      cplexPath=cplexPath, 
                                      numSolutions = numSolutions, 
                                      limitPop = limitPop, 
                                      poolIntensity = poolIntensity, 
                                      poolReplace = poolReplace){
      
      setwd(gsub(pattern = "cplex", replacement = "", x = cplexPath))
      data = commandFileName
      write(paste0("read ", inputFileName), data)
      write(paste0("set mip tolerances mipgap ",mipGap), data, append = TRUE)
      write(paste0("set mip pool relgap ",relGap), data, append = TRUE)
      write(paste0("set mip pool replace ", poolReplace), data, append = TRUE)
      write(paste0("set mip pool intensity ",poolIntensity), data, append = TRUE)
      write(paste0("set timelimit ",timelimit=timelimit), data, append = TRUE)
      write(paste0("set mip limits populate ", limitPop), data, append = TRUE)
      write(paste0("set mip pool capacity ", numSolutions), data, append = TRUE)
      write(paste0("populate"), data, append = TRUE)
      write(paste0("write ", outputFileName, " sol all"), data, append = TRUE)
      write("quit", data, append = TRUE)
      
    }
    
    # Writing and executing the cplex commands
    setwd(gsub(pattern = "cplex", replacement = "", x = cplexPath))
    writeCplexCommandFile(commandFileName = "cplexCommand.txt", inputFileName, outputFileName, mipGap=mipGap, relGap = relGap, timelimit=timelimit, cplexPath = cplexPath, numSolutions = numSolutions, limitPop = limitPop, poolIntensity = poolIntensity, poolReplace = poolReplace)
    if(file.exists(outputFileName)){
      file.remove(outputFileName)
    }
    system(paste0(cplexPath, " -f cplexCommand.txt"))
    
  }
  
  # Function that creates and runs the ILP problem
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
  
  # Function that reads the ILP solutions
  createILPBitstringAll<- function(cplexSolutionFileName,
                                   y_vector,
                                   binary_variables){
    
    #
    bitStringAll = list()
    
    #
    readOutBinariesAll<- function(cplexSolutionFileName, binary_variables){
      
      cplexSolutionData <- xmlParse(cplexSolutionFileName)
      cplexSolution <- xmlToList(cplexSolutionData, simplify = TRUE)
      binariesAll <- list()
      for(i in 1:length(cplexSolution)){
        
        if(all(names(cplexSolution[[i]])=="variable")){
          
          bitstring = c()
          binaries = c()
          for(j in 1:length(cplexSolution[[i]])){
            
            binaries[which(binary_variables[[1]]==cplexSolution[[i]][[j]][1])] = cplexSolution[[i]][[j]][3]
            
          }
          
          binariesAll[[length(binariesAll)+1]] <- binaries
          
        }
        
      }
      
      return(binariesAll)
      
    }
    
    #
    binaries = readOutBinariesAll(cplexSolutionFileName, binary_variables)
    for(i in 1:length(binaries)){
      bitString <- c()
      bitString = binaries[[i]][1:length(y_vector)]
      bitStringAll[[length(bitStringAll)+1]] <- bitString
    }
    
    #
    bitStringILPAll = bitStringAll
    for(i in 1:length(bitStringILPAll)){
      
      bitStringILPAll[[i]] <- round(as.numeric(bitStringILPAll[[i]]))
      
    }
    bitStringILPAll <- unique(bitStringILPAll)
    idx <- c()
    for(i in 1:length(bitStringILPAll)){
      
      if(sum(as.numeric(bitStringILPAll[[i]]))==0){
        
        idx <- c(idx, i)
        
      }
      
    }
    if(length(idx) > 0){
      
      bitStringILPAll <- bitStringILPAll[-idx]
      
    }
    
    return(bitStringILPAll)
    
  }
  
  # Initializing auxilliary objects for the analysis
  resILPAll <- list()
  exclusionList <- NULL
  cnolistReal <- cnolist
  tmpPath = getwd()
  
  # Checking cnolist class
  if ((class(cnolist)=="CNOlist")==TRUE){
    writeMIDAS(CNOlist = cnolist, filename = "tempMD.csv", timeIndices = c(1, 3), overwrite = TRUE)
    md = readMIDAS(MIDASfile = "tempMD.csv")
    file.remove("tempMD.csv")
    cnolist <- compatCNOlist(object = cnolist)
  }
  
  # Form of the objective function
  if(!(method%in%c("quadratic", "linear"))){
    print("'method' variable should be quadratic or linear. Setting method = 'qadratic'")
    method = "quadratic"
  }
  
  # Penalty factors as decimals
  options(scipen = 10)
  
  # Creating CNOlist object to use from the midas md table and checking if signals in the cnolist and pkn match
  cnolist = makeCNOlist(md,FALSE)
  
  startILP = Sys.time()
  print("Creating LP file and running ILP...")
  resILP = createAndRunILP(model, md, cnolistReal, accountForModelSize = TRUE, sizeFac = sizeFac, mipGap=mipGap, relGap=relGap, timelimit=timelimit, cplexPath = cplexPath, method = method, 
                           numSolutions = numSolutions, limitPop = limitPop, poolIntensity = poolIntensity, poolReplace = poolReplace)
  endILP = Sys.time()
  
  setwd(tmpPath)
  
  return(resILP)
  
}

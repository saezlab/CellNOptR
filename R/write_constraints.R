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

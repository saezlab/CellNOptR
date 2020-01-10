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
  
  ## constraint numbers refer to Mitsos et al., PloS Comp Biol. 2009
  ## in th .lp file, the constraint numbers are written for debugging
  
  constraintsVector <- writeConstraint3(model, midasExperimentPart)
  
  # constraints 4, 5 & 7
  reducedMidas <- midasTreatmentPart
  counter <- length(constraintsVector)+1
  for(i in 1:dim(reducedMidas)[1]){
    for(j in 1:length(model$reacID)){
      y_index <- j
      if(length(reaction_sets[[j]][[1]])!=0){ ## constraint 4
        for(k in 1:length(reaction_sets[[j]][[1]])){
          z_index <- length(model$reacID) + 
            length(model$namesSpecies)*dim(reducedMidas)[1] + 
            (i-1)*length(model$reacID) + j 
          x_index <- length(model$reacID) + 
            (i-1)*length(model$namesSpecies) + 
            which(model$namesSpecies==reaction_sets[[j]][[1]][k])
          newConstraint <- paste0("c", counter, ":\t", "xb",z_index, 
                                 " - xb",x_index, " <= 0 \t \t \\ c4")
          constraintsVector[counter] <- newConstraint
          counter= counter+1
        }
      }
      if(length(reaction_sets[[j]][[2]])!=0){ ## constraint 5
        for(k in 1:length(reaction_sets[[j]][[2]])){
          z_index <- length(model$reacID) + 
            length(model$namesSpecies)*dim(reducedMidas)[1] + 
            (i-1)*length(model$reacID) + j 
          x_index <- length(model$reacID) + (i-1)*length(model$namesSpecies) + 
            which(model$namesSpecies==reaction_sets[[j]][[2]][k])
          newConstraint <- paste0("c", counter, ":\t", "xb",z_index, " + xb",
                                 x_index, " <= 1\t \t \\ c5")
          constraintsVector[counter] <- newConstraint
          counter= counter+1
        }
      }
      if(length(reaction_sets[[j]][[3]])!=0){ ## constraint 7
        for(k in 1:length(reaction_sets[[j]][[3]])){
          speciesIsInhibitedInThisExperiment <- FALSE
          if(length(which(colnames(reducedMidas) ==
                          paste0("TR:",reaction_sets[[j]][[3]][k], "i")))==1){
            speciesIsInhibitedInThisExperiment <- 
              (reducedMidas[i,which(colnames(reducedMidas) ==
                                      paste0("TR:",reaction_sets[[j]][[3]][k], 
                                             "i"))] == 1)
          }
          
          if(!speciesIsInhibitedInThisExperiment){
            z_index <- length(model$reacID) + 
              length(model$namesSpecies)*dim(reducedMidas)[1] + 
              (i-1)*length(model$reacID) + j
            x_index <- length(model$reacID) + 
              (i-1)*length(model$namesSpecies) + 
              which(model$namesSpecies==reaction_sets[[j]][[3]][k])
            newConstraint <- paste0("c", counter, ":\t", "xb",z_index, " - xb",
                                   x_index, " <= 0 \t \t \\ c7")
            constraintsVector[counter] <- newConstraint
            counter= counter+1
          }
          
        }
      }
      
      
      #constraint 6
      z_index <- length(model$reacID) + 
        length(model$namesSpecies)*dim(reducedMidas)[1]+ 
        length(model$reacID)*(i-1) + j
      x_reactant_index <- c()
      x_inhibitor_index <- c()
      if(length(reaction_sets[[j]][[1]])!=0){## reactant set is not empty
        for(k in 1:length(reaction_sets[[j]][[1]])){
          x_reactant_index <- append(x_reactant_index, length(model$reacID) + 
                                      (i-1)*length(model$namesSpecies) +
                                      which(model$namesSpecies==
                                              reaction_sets[[j]][[1]][k]))
        }
      }
      
      if(length(reaction_sets[[j]][[2]])!=0){## inhibitors set is not empty
        for(k in 1:length(reaction_sets[[j]][[2]])){
          x_inhibitor_index  <- append(x_inhibitor_index, 
                                      length(model$reacID) + 
                                        (i-1)*length(model$namesSpecies) + 
                                        which(model$namesSpecies==
                                                reaction_sets[[j]][[2]][k]))
        }
      }
      
      newConstraint <- paste0("xb", y_index, " - ", "xb", z_index)
      reactants_sum <- paste0(" <= ", length(x_reactant_index))
      inhibitors_sum <- ""
      if(!is.null(x_reactant_index)){
        binary_sum <- ""
        for(l in 1:length(x_reactant_index)){
          binary_sum <- paste0(binary_sum, " + xb", x_reactant_index[l])
        }
        reactants_sum <- paste0( binary_sum, " <= " , length(x_reactant_index))
      }
      if(!is.null(x_inhibitor_index)){
        binary_sum <- ""
        for(l in 1:length(x_inhibitor_index)){
          binary_sum <- paste0(binary_sum, " - xb", x_inhibitor_index[l])
        }
        inhibitors_sum <- binary_sum
      }
      
      newConstraint <- paste0("c", counter,":\t", newConstraint, 
                             inhibitors_sum, reactants_sum, " \t \t \\ c6")
      constraintsVector[counter] <- newConstraint
      counter= counter+1
      ## --- End Constraint 6
    }
    
    #new constraint 8
    for(j in 1:length(model$namesSpecies)){
      flag_vector <- c()
      newConstraint <- paste0("xb", 
                             length(y_vector)+
                               (i-1)*length(model$namesSpecies)+j)
      for(k in 1:length(model$reacID)){
        
        if(reaction_sets[[k]][3] == model$namesSpecies[j]){
          flag_vector <- c(flag_vector,k)
          #}
        }
      }
      added_reactions <- ""
      if(!is.null(flag_vector)){
        for(k in 1:length(flag_vector)){
          added_reactions <- paste0(added_reactions,
                                   "- xb", 
                                   length(y_vector) + 
                                     (length(model$namesSpecies))*
                                     dim(reducedMidas)[1] + 
                                     (i-1)*length(model$reacID) + 
                                     flag_vector[k])
        }
        newConstraint <- paste0("c", counter, ": ", newConstraint, 
                               added_reactions, " <= 0" ," \t \t \\ c 8 new")
        constraintsVector[counter] <- newConstraint
        counter= counter+1  
      }
    }##end new 8
  }
  
  return(constraintsVector)
}
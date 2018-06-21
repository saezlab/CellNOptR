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
simulatorT0<-function(CNOlist,model,simList,indexList){

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }

    # This simulator is dedicated to the time0 case.
    # It creates a new structure called CNOlistT0, which is a subset of a CNOlist
    # It also sets the valueInhibitors and valueStimuli to zero.
    mi = CNOlist@inhibitors
    ms = CNOlist@stimuli
    # only inhibitors and stimuli are required but need to convert this into a
    # a valid CNOlist class, so we must provide all values

    CNOlistT0 = CNOlist
    # not needed, see the extension by AG below: 
    # CNOlistT0@inhibitors = matrix(1, dim(mi)[1], dim(mi)[2])
    # CNOlistT0@stimuli    = matrix(0, dim(ms)[1], dim(ms)[2])

    # !! hack. SimulatorT1 expect the inhibitors values to be 0 or 1.
    # Then, it flips them. So, the 0-values are flipped to 1. finally,
    # the values that are equal to 1 are set to NA...
    # Here, since this is T0, all inhibitors are set to 0. simulatorT1 will
    # therefore flip them to 1 and finally NA. So, we set all inhibitors to 1.
	
    # Need to be very careful if simulatorT1 changes
	
    # extension by AG 21.06.2018
    # As of 21.06.2018 the above proceedure is unneccesary. In simulatorT1.C, 
    # when mode == 0, all the inhibitors are initialised to 1 anyways. They are flipped
    # to 0 and therefore the simulated nodes will be innitialised to 0 value. But in the 
    # iterations of the simulation, (checks for mode == 0) the simulated nodes are
    # NOT overwritten by the inhibitory value, 
    # the value of the inhibitory nodes are calculated from the incoming edges.
    #
    # implementation of the permanentInhibitions / permanentStimuli nodes:
    # we change the mode 0 to mode 1 to pass stimuli and inhibitory values to the C simulator.
    # Non-permanent Stimuli are set to 0.
    # Non-permanent Inhibitors are set to 0, permanent inhibitors are set to 1.
    
    CNOlistT0 = CNOlist
    CNOlistT0@inhibitors = mi
    CNOlistT0@inhibitors[CNOlistT0@permanentInhibitors==0] = 0
    CNOlistT0@stimuli    = ms
    CNOlistT0@stimuli[CNOlistT0@permanentStimuli==0] = 0
    
    # Finally run the simulator with the particular set of experiments at t0
    newInput = simulatorT1(CNOlistT0, model, simList, indexList, mode=1)
    
    

    
    return(newInput)
    }


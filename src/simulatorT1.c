//============================================================================
// Name        : simulatorT1.c
// Author      : Aidan MacNamara
// Version     : 0.1
// Copyright   :
// Description :
//============================================================================

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>


SEXP simulatorT1 (

	SEXP nStimuli_in,
	SEXP nInhibitors_in,
	SEXP nCond_in,
	SEXP nReacs_in,
	SEXP nSpecies_in,
	SEXP nMaxInputs_in,
	
	SEXP finalCube_in,
	SEXP ixNeg_in,
	SEXP ignoreCube_in,
	SEXP maxIx_in,		
	
	SEXP indexSignals_in, 
	SEXP indexStimuli_in, 
	SEXP indexInhibitors_in,
			
	SEXP valueInhibitors_in,
	SEXP valueStimuli_in
) {
	
	SEXP simResults;
	
	int counter = 0;
	int i = 0;
	int j = 0;
	int curr_max = 0;
	int or_max = 0;
	int selection[40];
	int selCounter = 0;
	double *rans;
	
	int nStimuli = INTEGER(nStimuli_in)[0];
	int nInhibitors = INTEGER(nInhibitors_in)[0];
	int nCond = INTEGER(nCond_in)[0];
	int nReacs = INTEGER(nReacs_in)[0];
	int nSpecies = INTEGER(nSpecies_in)[0];
	int nMaxInputs = INTEGER(nMaxInputs_in)[0];

	counter = 0;
	int *maxIx;
	maxIx = (int*) malloc(nReacs * sizeof(int));
	for (i = 0; i < nReacs; i++) {
		maxIx[i] = INTEGER(maxIx_in)[counter++];
	}
	
	counter = 0;
	int *indexStimuli;
	indexStimuli = (int*) malloc(nStimuli * sizeof(int));
	for (i = 0; i < nStimuli; i++) {
		indexStimuli[i] = INTEGER(indexStimuli_in)[counter++];
	}

	counter = 0;
	int *indexInhibitors;
	indexInhibitors = (int*) malloc(nInhibitors * sizeof(int));
	for (i = 0; i < nInhibitors; i++) {
		indexInhibitors[i] = INTEGER(indexInhibitors_in)[counter++];
	}
	
	counter=0;
	int **finalCube;
	finalCube = (int**) malloc(nReacs * sizeof(int*));
	for (i = 0; i < nReacs; i++) {
	    finalCube[i] = (int*) malloc(nMaxInputs * sizeof(int));
	    for (j = 0; j < nMaxInputs; j++) {
			finalCube[i][j] = INTEGER(finalCube_in)[counter++];
	    }
	}
	
	counter=0;
	int **ixNeg;
	ixNeg = (int**) malloc(nReacs * sizeof(int*));
	for (i = 0; i < nReacs; i++) {
	    ixNeg[i] = (int*) malloc(nMaxInputs * sizeof(int));
	    for (j = 0; j < nMaxInputs; j++) {
			ixNeg[i][j] = INTEGER(ixNeg_in)[counter++];
	    }
	}
	
	counter=0;
	int **ignoreCube;
	ignoreCube = (int**) malloc(nReacs * sizeof(int*));
	for (i = 0; i < nReacs; i++) {
	    ignoreCube[i] = (int*) malloc(nMaxInputs * sizeof(int));
	    for (j = 0; j < nMaxInputs; j++) {
			ignoreCube[i][j] = INTEGER(ignoreCube_in)[counter++];
	    }
	}
	
	counter=0;
	int **valueInhibitors;
	valueInhibitors = (int**) malloc(nCond * sizeof(int*));
	for (i = 0; i < nCond; i++) {
	    valueInhibitors[i] = (int*) malloc(nInhibitors * sizeof(int));
	    for (j = 0; j < nInhibitors; j++) {
			valueInhibitors[i][j] = INTEGER(valueInhibitors_in)[counter++];
	    }
	}
	
	counter=0;
	int **valueStimuli;
	valueStimuli = (int**) malloc(nCond * sizeof(int*));
	for (i = 0; i < nCond; i++) {
	    valueStimuli[i] = (int*) malloc(nStimuli * sizeof(int));
	    for (j = 0; j < nStimuli; j++) {
			valueStimuli[i][j] = INTEGER(valueStimuli_in)[counter++];
	    }
	}
	
	
	//============================================================================
	
	// fill end_ix - how many reactions have each species as output
	int end_ix[nSpecies];
	int count_species=0;
	for(i = 0; i < nSpecies; i++) {
		for(j = 0; j < nReacs; j++) {
			if(i == maxIx[j]) {
				count_species++;
			}
		}
		end_ix[i] = count_species;
		count_species = 0;
	}
	
	// see stop conditions
	float test_val = 1e-3;
	
	// create an initial values matrix
	int init_values[nCond][nSpecies];
	for(i = 0; i < nCond; i++) {
		for(j = 0; j < nSpecies; j++) {
			init_values[i][j] = 2;
		}
	}
	
	// set the initial values of the stimuli
	for(i = 0; i < nCond; i++) {
		for(j = 0; j < nStimuli; j++) {
			init_values[i][indexStimuli[j]] = valueStimuli[i][j];
		}
	}
	
	// flip and redefine inhibitors
	if(nInhibitors) {
		for(i = 0; i < nCond; i++) {
			for(j = 0; j < nInhibitors; j++) {
				valueInhibitors[i][j] = 1 - valueInhibitors[i][j];
				if(valueInhibitors[i][j] == 1) {
					valueInhibitors[i][j] = 2;
				}
			}
		}
	}
	
	// set the initial values of the inhibitors
	for(i = 0; i < nCond; i++) {
		for(j = 0; j < nInhibitors; j++) {
			init_values[i][indexInhibitors[j]] = valueInhibitors[i][j];
		}
	}
	
	// initialize main loop
	int output_prev[nCond][nSpecies];
	int new_input[nCond][nSpecies];
	memcpy(new_input, init_values, sizeof(new_input));
	int term_check_1 = 1;
	float term_check_2 = 1;
	int count = 1;
	int diff;
	
	// define the temp data
	int temp_store[nCond * nReacs][nMaxInputs];
	
	
	//============================================================================
	
	// start simulation loop
	while(term_check_1 && term_check_2) {
		
		// copy to outputPrev
		memcpy(output_prev, new_input, sizeof(output_prev));
		
		// fill temp store
		// this is different to R version, through a single loop
		// with conditions
		int track_cond = 0; // track condition
		int track_reac = 0; // track reaction
		for(i = 0; i < nCond * nReacs; i++) {
			for(j = 0; j < nMaxInputs; j++) {
				// initial values of each input
				temp_store[i][j] = output_prev[track_cond][finalCube[track_reac][j]];

				if(ignoreCube[track_reac][j]) {
					temp_store[i][j] = 2;
				}
				if(ixNeg[track_reac][j]) {
					// flip the values of the neg inputs
					if(temp_store[i][j] == 0) {temp_store[i][j] = 1;}
					else if(temp_store[i][j] == 1) {temp_store[i][j] = 0;}
				}
			//		Rprintf("%d\n", temp_store[i][j]);

			}
			
			track_cond++;
			if((track_cond == nCond)) {
				track_cond = 0;
				track_reac++;
			}
		}
	
		// compute the AND gates (find the min 0/1 of each row)
		int output_cube[nCond][nReacs]; // declare output_cube
		
		int current_min;
		int dial_reac = 0;
		int dial_cond = 0;
		for(i = 0; i < nCond * nReacs; i++) {
			current_min = temp_store[i][0];
			for(j = 1; j < nMaxInputs; j++) {
				if(temp_store[i][j] < current_min) {current_min = temp_store[i][j];}
			}
			output_cube[dial_cond][dial_reac] = current_min;
			dial_cond++;
			if(dial_cond==nCond) {dial_cond = 0; dial_reac++;}
		}
		
		// compute the OR gates and reinitialize new_input
		for(i = 0; i < nCond; i++) {
			for(j = 0; j < nSpecies; j++) {
				new_input[i][j] = 2;
			}
		}
		
		// declare vector to store 'selection' (R)
		selCounter = 0;
		for(int s = 0; s < nSpecies; s++) {
			// is the species an output for any reactions?
			if(end_ix[s]) {

				// find reactions with this species as output
				// add equivalent output_cube data to new_input
				for(int a = 0; a < nReacs; a++) {
					if(s == maxIx[a]) {selection[selCounter] = a; selCounter++;}
				}
				// if the species is an output for a single reaction
				// it's a 1-1 mapping to new_input
				if(selCounter == 1) {
					for(int b = 0; b < nCond; b++) {
						new_input[b][s] = output_cube[b][selection[selCounter-1]];
					}
					selCounter = 0;
				}
				// else if species is output for > 1
				if(selCounter > 1) {
					for(i=0; i < nCond; i++) {
						or_max = 2;
						curr_max = 0;
						for(int p=0; p < selCounter; p++) {
							if(output_cube[i][selection[p]] >= curr_max && output_cube[i][selection[p]] < 2) {
								or_max = output_cube[i][selection[p]];
								curr_max = output_cube[i][selection[p]];
							}
						}
						new_input[i][s] = or_max;
					}
					selCounter = 0;
				}
			}
		}
		
		// reset the stimuli
		for(i = 0; i < nCond; i++) {
			for(j = 0; j < nStimuli; j++) {
				curr_max = valueStimuli[i][j];
				if(new_input[i][indexStimuli[j]] > curr_max && new_input[i][indexStimuli[j]] < 2) {
					curr_max = new_input[i][indexStimuli[j]];
				}
				new_input[i][indexStimuli[j]] = curr_max;
			}
		}
		
		// reset the inhibitors
		for(i = 0; i < nCond; i++) {
			for(j = 0; j < nInhibitors; j++) {
				if(valueInhibitors[i][j] == 0) {
					new_input[i][indexInhibitors[j]] = 0;
				}
			}
		}
		
		// set 'NAs' (2s) to 0
		for(i = 0; i < nCond; i++) {
			for(j = 0; j < nSpecies; j++) {
				if(new_input[i][j] == 2) {new_input[i][j] = 0;}
				if(output_prev[i][j] == 2) {output_prev[i][j] = 0;}
			}
		}
		

        term_check_1 = 0;
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nSpecies; j++) {
                diff = abs((new_input[i][j] - output_prev[i][j]));
                if (diff > test_val){
                    term_check_1 = 1;
                    break;  /*  no need to keep going checking other values if 
                                one is greater than test_val */
                }
            }
        }
        /*term_check_1 = !(abs(diff) < test_val);*/
        term_check_2 = (count < (nSpecies * 1.2));
		count++;		
		
	}
	
	// set non-resolved bits to 2 (NA)
	for(i = 0; i < nCond; i++) {
		for(j = 0; j < nSpecies; j++) {
			if(new_input[i][j] != output_prev[i][j])
				new_input[i][j] = 2;
		}
	}
	
	PROTECT(simResults = allocMatrix(REALSXP, nCond, nSpecies));
	rans = REAL(simResults);
	for(i = 0; i < nCond; i++) {
		for(j = 0; j < nSpecies; j++) {
			if(new_input[i][j] == 2) rans[i + nCond*j] = NA_REAL;
			else rans[i + nCond*j] = new_input[i][j];
		}
	} 

	
	free(maxIx);
	free(indexStimuli);
	free(indexInhibitors);

	for (i = 0; i < nReacs; i++) {
		free(finalCube[i]);
	}
	free(finalCube);
	
	for (i = 0; i < nReacs; i++) {
		free(ixNeg[i]);
	}
	free(ixNeg);
	
	for (i = 0; i < nReacs; i++) {
		free(ignoreCube[i]);
	}
	free(ignoreCube);
	
	for (i = 0; i < nCond; i++) {
		free(valueInhibitors[i]);
	}
	free(valueInhibitors);
	
	for (i = 0; i < nCond; i++) {
		free(valueStimuli[i]);
	}
	free(valueStimuli);

	UNPROTECT(1);
	return simResults;

}

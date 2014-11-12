/* ============================================================================
 *  Name        : simulatorTN.c
 *  Author      : Thomas Cokelaer based on simulatorT1 (Aidan MacNamara)
 *  Version     : 0.1
 *  Copyright   :
 *  Description : Requires slightly more inputs (e.g., previous results at TN-1,
 *   interMat that may be
 *   replaced by the contents of other variables, and times variables.). The
 *   main difference being in the requirement of some initialisation. This
 *   implementation is similar to the R version to ease debugging. Once stable,
 *   it may diverge from the R version around 1800).
 *
 * =========================================================================== */
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>


/* keep this NA > 1 and integer */
#define NA 100

SEXP simulatorTN (

    SEXP nStimuli_in,
    SEXP nInhibitors_in,
    SEXP nCond_in,
    SEXP nReacs_in,
    SEXP nSpecies_in,
    SEXP nSignals_in,
    SEXP nMaxInputs_in,
    SEXP nTimes_in,


    SEXP timeIndex_in,
    SEXP times_in,
    SEXP interMat_in,

    SEXP prevSimResults_in,

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
    int k = 0;
    int s = 0;
    int curr_max = 0;
    int or_max = 0;
    int selection[40];
    int selCounter = 0;
    double *rans;


    /* variables */
    int nStimuli = INTEGER(nStimuli_in)[0];
    int nInhibitors = INTEGER(nInhibitors_in)[0];
    int nCond = INTEGER(nCond_in)[0];
    int nReacs = INTEGER(nReacs_in)[0];
    int nSpecies = INTEGER(nSpecies_in)[0];
    int nSignals = INTEGER(nSignals_in)[0];
    int nMaxInputs = INTEGER(nMaxInputs_in)[0];
    int nTimes = INTEGER(nTimes_in)[0];

    int timeIndex = INTEGER(timeIndex_in)[0];

    /* variables to read the inputs*/
    int *times;
    int *maxIx;
    int *indexStimuli;
    int *indexInhibitors;
    int *indexSignals;
    int **interMat;
    int **prevSimResults;
    int **finalCube;
    int **ixNeg;
    int **ignoreCube;
    int **valueInhibitors;
    int **valueStimuli;

    int *end_ix;
    int count_species=0;

    /* see stop conditions */
    float test_val = 1e-3;
    int **output_prev;
    int **new_input;
    int **first_iter;
/*    int **output_cube;*/
  /* int **temp_store;*/

    int track_cond = 0; /* track condition */
    int track_reac = 0; /* track reaction */

    int current_min;
    int dial_reac = 0;
    int dial_cond = 0;

    int nreacsTN = 0;
    int *reacsTN;

    int outNode, r;
    int a,b,p;
    int term_check_1 = 1;
    float term_check_2 = 1;
    int count = 1;
    int diff;

    counter = 0;
    times = (int*) malloc(nTimes * sizeof(int));
    for (i = 0; i < nTimes; i++) {
        times[i] = INTEGER(times_in)[counter++];
    }

    counter = 0;
    maxIx = (int*) malloc(nReacs * sizeof(int));
    for (i = 0; i < nReacs; i++) {
        maxIx[i] = INTEGER(maxIx_in)[counter++];
    }

    counter = 0;
    indexStimuli = (int*) malloc(nStimuli * sizeof(int));
    for (i = 0; i < nStimuli; i++) {
        indexStimuli[i] = INTEGER(indexStimuli_in)[counter++];
    }

    counter = 0;
    indexInhibitors = (int*) malloc(nInhibitors * sizeof(int));
    for (i = 0; i < nInhibitors; i++) {
        indexInhibitors[i] = INTEGER(indexInhibitors_in)[counter++];
    }

    counter = 0;
    indexSignals = (int*) malloc(nSignals * sizeof(int));
    for (i = 0; i < nSignals; i++) {
        indexSignals[i] = INTEGER(indexSignals_in)[counter++] ; 
    }

    counter=0;
    interMat = (int**) malloc(nSpecies * sizeof(int*));
    for (i = 0; i < nSpecies; i++) {
        interMat[i] = (int*) malloc(nReacs * sizeof(int));
        for (j = 0; j < nReacs; j++) {
            interMat[i][j] = INTEGER(interMat_in)[counter++];
        }
    }

    /* previous SimResults data */
    counter=0;
    prevSimResults = (int**) malloc(nCond * sizeof(int*));
    for (i = 0; i < nCond; i++) {
        prevSimResults[i] = (int*) malloc(nSpecies * sizeof(int));
        for (j = 0; j < nSpecies; j++) {
            prevSimResults[i][j] = INTEGER(prevSimResults_in)[counter++];
        }
    }

    counter=0;
    finalCube = (int**) malloc(nReacs * sizeof(int*));
    for (i = 0; i < nReacs; i++) {
        finalCube[i] = (int*) malloc(nMaxInputs * sizeof(int));
        for (j = 0; j < nMaxInputs; j++) {
            finalCube[i][j] = INTEGER(finalCube_in)[counter++];
        }
    }

    counter=0;
    ixNeg = (int**) malloc(nReacs * sizeof(int*));
    for (i = 0; i < nReacs; i++) {
        ixNeg[i] = (int*) malloc(nMaxInputs * sizeof(int));
        for (j = 0; j < nMaxInputs; j++) {
            ixNeg[i][j] = INTEGER(ixNeg_in)[counter++];
        }
    }

    counter=0;
    ignoreCube = (int**) malloc(nReacs * sizeof(int*));
    for (i = 0; i < nReacs; i++) {
        ignoreCube[i] = (int*) malloc(nMaxInputs * sizeof(int));
        for (j = 0; j < nMaxInputs; j++) {
            ignoreCube[i][j] = INTEGER(ignoreCube_in)[counter++];
        }
    }

    counter=0;
    valueInhibitors = (int**) malloc(nCond * sizeof(int*));
    for (i = 0; i < nCond; i++) {
        valueInhibitors[i] = (int*) malloc(nInhibitors * sizeof(int));
        for (j = 0; j < nInhibitors; j++) {
            valueInhibitors[i][j] = INTEGER(valueInhibitors_in)[nCond*j+i];
        }
    }

    counter=0;
    valueStimuli = (int**) malloc(nCond * sizeof(int*));
    for (i = 0; i < nCond; i++) {
        valueStimuli[i] = (int*) malloc(nStimuli * sizeof(int));
        for (j = 0; j < nStimuli; j++) {
            valueStimuli[i][j] = INTEGER(valueStimuli_in)[nCond*j+i];
        }
    }

    /* no need to set the initial values for those variables*/
  /*  temp_store = (int**) malloc(nCond*nReacs * sizeof(int*));


    for (i = 0; i < nCond*nReacs; i++) {
        temp_store[i] = (int*) malloc(nMaxInputs * sizeof(int));
    }*/

    output_prev = (int**) malloc(nCond*sizeof(int*));
    new_input = (int**) malloc(nCond*sizeof(int*));
    first_iter = (int**) malloc(nCond*sizeof(int*));
/*    output_cube = (int**) malloc(nCond*sizeof(int*));*/
    for (i = 0; i < nCond; i++) {
        output_prev[i] = (int*) malloc(nSpecies * sizeof(int*));
        new_input[i] = (int*) malloc(nSpecies * sizeof(int*));
        first_iter[i] = (int*) malloc(nSpecies * sizeof(int*));
         /*note that here we need nReacs instead of nSpecies*/
/*        output_cube[i] = (int*) malloc(nReacs * sizeof(int*));*/
    }


    /*============================================================================*/

    /* fill end_ix - how many reactions have each species as output*/
    end_ix = (int*) malloc(nSpecies * sizeof(int));
    count_species=0;
    for(i = 0; i < nSpecies; i++) {
        for(j = 0; j < nReacs; j++) {
            if(i == maxIx[j]) {
                count_species++;
            }
        }
        end_ix[i] = count_species;
        count_species = 0;
    }

    /* see stop conditions */
    test_val = 1e-3;


    /* First iteration before main loop. Different from T1 !*/
    int temp_store[nCond * nReacs][nMaxInputs];
  int output_cube[nCond][nReacs]; 

    /* get back previous results */
    for (i=0; i<nCond ;i++){
        for (j=0; j<nSpecies ;j++){
            new_input[i][j] = prevSimResults[i][j];
         }
     }

   // Rprintf("AAA\n");


    for (i=0; i<nCond ;i++){
        for (j=0; j<nSpecies ;j++){
            output_prev[i][j] = new_input[i][j];
         }
     }


    /* 1. need to compute tempStore */

    /* R code
      tempStore<-apply(simList$finalCube,2,function(x){return(outputPrev[,x])})
      tempIxNeg<-apply(simList$ixNeg,2,filltempCube)
      tempIgnore<-apply(simList$ignoreCube,2,filltempCube)
      tempStore[tempIgnore]<-NA
      tempStore[tempIxNeg]<-1-tempStore[tempIxNeg]
    */
    track_cond = 0; 
    track_reac = 0; 
    for(i = 0; i < nCond * nReacs; i++) {
        for(j = 0; j < nMaxInputs; j++) {
            /* initial values of each input */
            temp_store[i][j] = output_prev[track_cond][finalCube[track_reac][j]];

            if(ignoreCube[track_reac][j]) {
                temp_store[i][j] = 2;
            }
            if(ixNeg[track_reac][j]) {
                /* flip the values of the neg inputs */
                if(temp_store[i][j] == 0) {temp_store[i][j] = 1;}
                else if(temp_store[i][j] == 1) {temp_store[i][j] = 0;}
            }

        }
        track_cond++;
        if((track_cond == nCond)) {
            track_cond = 0;
            track_reac++;
        }
    }
   
    /* 2. Need to compute OutputCube
     R code :
       outputCube <- apply(tempStore, 1, minNA)
       outputCube<-matrix(outputCube, nrow=nCond,ncol=nReacs)
    */


    dial_reac = 0;
    dial_cond = 0;
    for(i = 0; i < nCond * nReacs; i++) {
        current_min = temp_store[i][0];
        for(j = 1; j < nMaxInputs; j++) {
		
			/* if statement below is for AND gates with any NA (2) input
			 in this case output should always be NA*/
			if(temp_store[i][j] == 2 && ignoreCube[dial_reac][j] == 0) {
				current_min = 2;
				break;
			}
			else if(temp_store[i][j] < current_min) {current_min = temp_store[i][j];}
		
		}
        output_cube[dial_cond][dial_reac] = current_min;
        dial_cond++;
        if(dial_cond==nCond) {dial_cond = 0; dial_reac++;}
    }


    /* 3. figure out the reacsTN */
    counter = 0;
    for (i=0; i<nTimes; i++){

        if (timeIndex - 1 == times[i]){
            counter++;
        }
    }
    nreacsTN = counter;
    reacsTN = (int*) malloc(nreacsTN * sizeof(int));
    counter = 0 ;
    for (i=0; i<nTimes; i++){
        if (timeIndex - 1 == times[i]){
            reacsTN[counter] = i;
            counter++;
        }
    }


    /* scan the interMat looking for specific reactions (reacsTN) that have
     interMat positive (i.e, search for B in A=B) */
    for (i=0; i<nreacsTN; i++){
        for (j=0; j<nSpecies; j++){
            if (interMat[j][reacsTN[i]]>0){
                /* outNode is the B in A=B */
                outNode = j;
                /* now we scan again the interMAt looking for A (multiple
                 solution possible) */
                for (r=0; r<nReacs; r++){
                    if (interMat[outNode][r]>0){
                        for (k=0; k<nCond; k++){
                            output_cube[k][r] = output_cube[k][reacsTN[i]];
                         }
                     }
                 }
             }
        }
    }

    /* 4. compute or gate 
     Note that we populate the matrix by columns and then rows which is not optimal... */
    selCounter = 0;
    for(s = 0; s < nSpecies; s++) {

        if(end_ix[s]) {

            /* find reactions with this species as output
             add equivalent output_cube data to new_input*/
            for(a = 0; a < nReacs; a++) {
                if(s == maxIx[a]) {selection[selCounter] = a; selCounter++;}
            }
            /* if the species is an output for a single reaction
             it's a 1-1 mapping to new_input */
            if(selCounter == 1) {
                for(b = 0; b < nCond; b++) {
                    new_input[b][s] = output_cube[b][selection[selCounter-1]];
                }
                selCounter = 0;
            }
            /* else if species is output for > 1 */
            if(selCounter > 1) {
                for(i=0; i < nCond; i++) {
                    or_max = 2;
                    curr_max = 0;
                    for(p=0; p < selCounter; p++) {
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
    /* 5. compute and gate */

    /* reset the stimuli */
    for(i = 0; i < nCond; i++) {
         for(j = 0; j < nStimuli; j++) {
             curr_max = valueStimuli[i][j];
             if(new_input[i][indexStimuli[j]] > curr_max && new_input[i][indexStimuli[j]] < 2) {
                 curr_max = new_input[i][indexStimuli[j]];
             }
             new_input[i][indexStimuli[j]] = curr_max;
         }
     }

    /*indexInhibitors = (int*) malloc(nInhibitors * sizeof(int));*/
    for (i=0; i<nCond; i++){
        for (j=0; j<nInhibitors; j++){
                    valueInhibitors[i][j] = 1 - valueInhibitors[i][j] ;
        }
    }
    for (i=0; i<nInhibitors; i++){
        for (j=0; j<nCond; j++){
            new_input[j][indexInhibitors[i]] = valueInhibitors[j][i] * new_input[j][indexInhibitors[i]] ;
        }
    }
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nSpecies; j++) {
            if (new_input[i][j]==2){
                new_input[i][j] = 0;
            }
        }
    }

    for (i=0; i<nCond ;i++){
        for (j=0; j<nSpecies ;j++){
            first_iter[i][j] = new_input[i][j];
         }
     }

    /*============================================================================*/


    /* reset output_cube */
    for (i=0; i<nCond; i++){
        for (j=0; j<nReacs; j++){
        output_cube[nCond][nReacs] = 0.;
        }
    }


    term_check_1 = 1;
    term_check_2 = 1;
    count = 1;


    /* start simulation loop*/
    while(term_check_1 && term_check_2) {

        /* copy to outputPrev*/
        for (i=0; i<nCond; i++){
            for(j=0; j<nSpecies; j++){
                output_prev[i][j] = new_input[i][j];
            }
        }


        /*
         fill temp store
         this is different to R version, through a single loop
         with conditions*/
        track_cond = 0; /* track condition*/
        track_reac = 0; /* track reaction*/
        for(i = 0; i < nCond * nReacs; i++) {
            for(j = 0; j < nMaxInputs; j++) {
                /* initial values of each input*/
                temp_store[i][j] = output_prev[track_cond][finalCube[track_reac][j]];

                if(ignoreCube[track_reac][j]) {
                    temp_store[i][j] = NA;
                }
                if(ixNeg[track_reac][j]) {
                    /* flip the values of the neg inputs*/
                    if(temp_store[i][j] == 0) {temp_store[i][j] = 1;}
                    else if(temp_store[i][j] == 1) {temp_store[i][j] = 0;}
                }

            }

            track_cond++;
            if((track_cond == nCond)) {
                track_cond = 0;
                track_reac++;
            }
        }


        /* compute the AND gates (find the min 0/1 of each row)*/

        dial_reac = 0;
        dial_cond = 0;
        for(i = 0; i < nCond * nReacs; i++) {
            current_min = temp_store[i][0];
            for(j = 1; j < nMaxInputs; j++) {
                if(temp_store[i][j] < current_min) {current_min = temp_store[i][j];}
            }

            output_cube[dial_cond][dial_reac] = current_min;
            dial_cond++;
            if(dial_cond==nCond) {dial_cond = 0; dial_reac++;}
        }

        /* compute the OR gates and reinitialize new_input*/
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nSpecies; j++) {
                new_input[i][j] = NA;
            }
        }

        /* declare vector to store 'selection' (R)*/
        selCounter = 0;
        for(s = 0; s < nSpecies; s++) {
            /* is the species an output for any reactions?*/
            if(end_ix[s]) {

                /* find reactions with this species as output
                 add equivalent output_cube data to new_input */
                for(a = 0; a < nReacs; a++) {
                    if(s == maxIx[a]) {selection[selCounter] = a; selCounter++;}
                }
                /* if the species is an output for a single reaction
                 it's a 1-1 mapping to new_input */
                if(selCounter == 1) {
                    for(b = 0; b < nCond; b++) {
                        new_input[b][s] = output_cube[b][selection[selCounter-1]];
                    }
                    selCounter = 0;
                }
                /* else if species is output for > 1 */
                if(selCounter > 1) {
                    for(i=0; i < nCond; i++) {
                        or_max = NA;
                        curr_max = 0;
                        for(p=0; p < selCounter; p++) {
                            if(output_cube[i][selection[p]] >= curr_max && output_cube[i][selection[p]] < NA) {
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

        /* reset the stimuli */
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nStimuli; j++) {
                curr_max = valueStimuli[i][j];
                if(new_input[i][indexStimuli[j]] > curr_max && new_input[i][indexStimuli[j]] < NA) {
                    curr_max = new_input[i][indexStimuli[j]];
                }
                new_input[i][indexStimuli[j]] = curr_max;
            }
        }

        /* reset the inhibitors */
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nInhibitors; j++) {
                if(valueInhibitors[i][j] == 0) {
                    new_input[i][indexInhibitors[j]] = 0;
                }
            }
        }

       for (i=0; i<nSpecies; i++){
           for (j=0; j<nreacsTN; j++){
              if (interMat[i][reacsTN[j]]>0){
                for (k=0; k<nCond; k++){
                  new_input[k][i] = first_iter[k][i];
                 }
              }
         }
      }

        /* set 'NAs' (2s) to 0 */
        for(i = 0; i < nCond; i++) {
            for(j = 0; j < nSpecies; j++) {
                if(new_input[i][j] == NA) {new_input[i][j] = 0;}
                if(output_prev[i][j] == NA) {output_prev[i][j] = 0;}
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

    } /* end of main loop/
     set non-resolved bits to 2 (NA)*/
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nSpecies; j++) {
            if(new_input[i][j] != output_prev[i][j])
                new_input[i][j] = NA;
        }
    }

     PROTECT(simResults = allocMatrix(REALSXP, nCond, nSpecies));
    rans = REAL(simResults);
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nSpecies; j++) {
            if(new_input[i][j] == NA) rans[i + nCond*j] = NA_REAL;
            else rans[i + nCond*j] = new_input[i][j];
        }
    }

/*
     PROTECT(simResults = allocMatrix(REALSXP, nCond, nSignals));
    rans = REAL(simResults);
    for(i = 0; i < nCond; i++) {
        for(j = 0; j < nSignals; j++) {
            if(new_input[i][indexSignals[j]] == NA) rans[i + nCond*j] = NA_REAL;
            else rans[i + nCond*j] = new_input[i][indexSignals[j]];
        }
    }

*/


    free(maxIx);
    free(indexStimuli);
    free(indexInhibitors);
    free(indexSignals);
    free(end_ix);
    free(times);
    free(reacsTN);

    for (i = 0; i < nSpecies; i++) {
        free(interMat[i]);
    }
    free(interMat);

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

    for (i = 0; i < nCond; i++) {
        free(prevSimResults[i]);
    }
    free(prevSimResults);


    for (i=0; i<nCond; i++){
        free(first_iter[i]);
        free(new_input[i]); 
       free(output_prev[i]);
    }
    free(first_iter);
    free(new_input);
    free(output_prev);




    UNPROTECT(1);
    return simResults;

}

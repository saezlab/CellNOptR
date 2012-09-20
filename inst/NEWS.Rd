\title{NEWS}
\name{CellNOptR NEWS}

\section{CHANGES IN VERSION 1.3.33: (2012-09-20)}{
    \itemize{
        \item{BUG FIX}
            \itemize{
                \item{writeMIDAS uses the new CNOlist class as an input. For
back compatibility, the input can also be the data frame output by makeCNOlist;
this option was buggy, which is now fixed.}
            }
}


\section{CHANGES IN VERSION 1.3.32: (2012-09-10)}{
    \itemize{
        \item{BUG FIX}
            \itemize{
                \item{plotModel can handle CNOlist=NULL as it used to be the case}
            }
        \item{CHANGES}
            \itemize{
                \item{simulatorT1 and TN have additional check for AND gates with any NA inputs.}
            }
    }
}

\section{CHANGES IN VERSION 1.3.31: (2012-09-06)}{
    \itemize{
        \item{BUG FIXES}
            \itemize{
                \item{typo in checkSignals in warning message.}
                \item{fixing bug (typo) in cutAndPlotResultsT1}
            }
   }
}

\section{CHANGES IN VERSION 1.3 to 1.3.30: (2012-08-30)}{
    \itemize{
        \item{BUG FIXES}
           \itemize{
                \item{fixed bug in writeDot function when only 1 stimuli present.}
                \item{cutNONC function recompte reacID on its own so that if compresion and
                  expansion are not performed, reacID are correct}
            }
        \item{CHANGES AND UI SIMPLIFICATION}
         \itemize{
            \item{major update that changes User Interface to uniformise the code for
                  users and developers (e.g., readSif function is now called readSIF and prep4Sim is prep4sim).
                  Similarly many arguments are not required anymore. For instance, indexList and simList are not required
                  anymore (for back compatibilty these 2 arguments are usually kept and
                    set to NULL).}
            \item{improvment of  plotModel function to have more graphviz options}
            \item{gaBinary default initBstring argument is set to ones.}
            \item{MAJOR refactoring to merge T1, T2, TN functions (cutAnPlot, gaBinary,  simulator).}
            \item{simulateTN  combines simulateT1 and simulateTN so only one function is required.}
            \item{gaBinaryT2 is deprecated. Use TN version instead.}
            \item{gaBinaryT1 and TN do no need the indexList and simList parameters anymore}
             \item{preprocessing now returns only the model. indices are not required anymore.
                      preprocessing can now handle model without cnolist provided (e.g. for expansion only)}
            \item{add a cutAndPlot dispatcher that calls the proper cutandPlotResultsTX function.}
             \item{cutSimList argument called bitString renamed bString (developers)}
            \item{computeScoreT1, cutSimList had an argument called bitString. it is now called
                   bString and this convention is used throughout the code (e.g., bStringT1 not
                   bitStringT1) simlist and indexlist are now optional are after the bString
                   argument.}
            \item{ show argument removed from cutAndPlot function family}
         }

        \item{NEW FEATURES}
            \itemize{
                \item{c Simulator implemented making the optimisation much faster.}
                \item{steady states extended from T1 and T2 to any time TN}
                \item{add a CNOlist class that is now used throughout the code (has a plot method that replaces plotCNOlist calls)}
                \item{add plotParams parameter in cutAndPlot family.}
                \item{hash table in gaBinaryT1 (maxSizeHashTable as optional argument)}
                \item{add cutModel function to make code more modular (for developers)}
                \item{add mapBack function}
                \item{add computeScore functions to make the code modular.}
                \item{conversion functions sif2graph, model2sif, graph2sif added to the package}
                \item{add writeMIDAS function}
                \item{add writeSIF function}
                \item{add prior knowledge bitstring in the GA algorithm (T1 and T2)}
                \item{add a test suite in tests directory}
            }
       }
    
}







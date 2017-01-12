---
layout: default
title: Home
---


# Welcome to the CellNOpt Documentation Page!

## Overview
**CellNOpt** (from CellNetOptimizer; a.k.a. CNO) is a software used for creating logic-based models of signal transduction networks using different logic formalisms (Boolean, Fuzzy, or differential equations). CellNOpt uses information on signaling pathways encoded as a Prior Knowledge Network, and trains it against high-throughput biochemical data to create cell-specific models.

CellNOpt is freely available under GPL license in R and Matlab languages. It can be also accessed through a python wrapper, and a Cytoscape plugin called [CytoCopter](http://www.cellnopt.org/cytocopter/index.html) provides a graphical user interface.

<img src="/cellnopt/public/index1.png" alt="Example result">


## CellNOpt Implementations
### _CellNOptR (R packages)_
A series of packages are available in R. The core CellNOpt is available on BioConductor web site: CellNOptR, revision 1.4.0. Newest and oldest version are also available in our [Downloads](http://www.ebi.ac.uk/saezrodriguez/cno/downloads.html) page.

[CellNOptR](http://www.bioconductor.org/packages/release/bioc/html/CellNOptR.html) contains the core functions as well as the boolean and steady states version. It implements the workflow described in Saez-Rodriguez et al Mol Sys Bio 2009, with extended capabilities for multiple time points.

[CNORdt](http://www.bioconductor.org/packages/release/bioc/html/CNORdt.html) is an extension that allows to train a Boolean model agains time-courses of data.

[CNORfuzzy](http://www.bioconductor.org/packages/release/bioc/html/CNORfuzzy.html) is an extension to CellNOptR that allows to handle continous values, using constrained fuzzy logic, as described in Morris et al Plos Comp Bio 2011.

[CNORode](http://www.bioconductor.org/packages/release/bioc/html/CNORode.html) is an ODE add-on to CellNOptR. It is based on the method of (Wittmann et al BMC Sys Bio 2009), also implemented in the tool Odefy (Krusiek et al BMC Bioinf 2010).

[CNORfeeder](http://www.bioconductor.org/packages/release/bioc/html/CNORfeeder.html) is an add-on to CellNOptR that permits to extend a network derived from literature with links derived in a strictly data-driven way and supported by protein-protein interactions as described in (Eduati et al Bioinformatics 2012).

### _MATLAB_
Some features of CellNOpt are also available as a MATLAB toolbox, along with the toolbox Q2LM to analyze models, [here](http://www.ebi.ac.uk/saezrodriguez/cno/matlab)

### _Python_
A Python package called cellnopt.wrapper provides a python interace to the R packages (CellNOptR, CNORode and CNORfuzzy). It uses rpy2 and is available on [Pypi](http://pypi.python.org/pypi/cellnopt.wrapper/). For more details see its [cellnopt.wrapper](http://www.ebi.ac.uk/~cokelaer/cellnopt/wrapper) page. In addition a pure Python version is developed on [github](http://github.com/cellnopt/cellnopt).

### _Cytoscape Plugin (CytoCopter)_
CytoCopteR is a Graphical User Interface designed as a Cytoscape plugin. It provides an interface to CellNOptR using Rserve. More information is available on [CytoCopter](http://www.cellnopt.org/cytocopter/index.html) page.


## Complementary Tools
### _MEIGO_
[MEIGO](http://www.iim.csic.es/~gingproc/meigo.html), a global optimization toolbox that includes a number of metaheuristic methods as well as a Bayesian inference method for parameter estimation, that can be applied to model training in CellNOpt. Available in R, Matlab, and Python. Presented in Egea et al BMC Bioinformatics, 214.

### _Caspo_
[Caspo](http://bioasp.github.io/caspo/) a Python toolbox based on Answer Set Programming to exactly and exhaustively train Boolean models as defined in CellNOpt’s Boolean steady state case. Presented in Guziolowski et al Bioinformatics, 2013 [link](http://bioinformatics.oxfordjournals.org/content/29/18/2320.long)

### _CoLoMoTo_
The [ColoMoTo](http://www.colomoto.org/) consortium nvolves other groups developing tools and methods for logic modelling. We have also jointly develop the standard SBML-qual (Chaouiya et al, BMC Syst Bio 2013) ([link](http://www.colomoto.org/)) that allows to exchange models within tools.

### _PHONEMeS_
[PHONEMeS](http://saezlab.github.io/PHONEMeS/) Toolbox dedicated to mass spectrometry analysis.


## References
+ Main reference describing CellNOpt, which can be used to cite it: C Terfve, T Cokelaer, A MacNamara, D Henriques, E Goncalves, MK Morris, M van Iersel, DA Lauffenburger, J Saez-Rodriguez. [CellNOptR: a flexible toolkit to train protein signaling networks to data using multiple logic formalisms](http://www.biomedcentral.com/1752-0509/6/133/abstract). BMC Systems Biology, 2012, 6:133 [PDF](http://www.biomedcentral.com/content/pdf/1752-0509-6-133.pdf)

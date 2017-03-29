---
layout: default
title: Home
---


# Welcome to the CellNOpt Documentation Page!

## Overview
**CellNOpt** (from CellNetOptimizer; a.k.a. CNO) is a software used for creating logic-based models of signal transduction networks using different logic formalisms (Boolean, Fuzzy, or differential equations). CellNOpt uses information on signaling pathways encoded as a Prior Knowledge Network, and trains it against high-throughput biochemical data to create cell-specific models.

CellNOpt is freely available under GPL license in R and Matlab languages. It can be also accessed through a python wrapper, and a Cytoscape plugin called [CytoCopter](http://www.cellnopt.org/cytocopter/index.html) provides a graphical user interface.

<img src="/cellnopt/public/index1.png" alt="Example result">

| **CellNOpt** is described in details in the following paper (more literature related to CellNOpt is available in the Publications     sections - below this page). Please use this reference to cite CellNOpt:  |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **C Terfve, T Cokelaer, A MacNamara, D Henriques, E Goncalves, MK Morris, M van Iersel, DA Lauffenburger, J Saez-Rodriguez.** [CellNOptR: a flexible toolkit to train protein signaling networks to data using multiple logic formalisms](http://www.biomedcentral.com/1752-0509/6/133/abstract). _BMC Systems Biology, 2012, 6:133_ [PDF](http://www.biomedcentral.com/content/pdf/1752-0509-6-133.pdf)  |

| We have also developed PHONEMeS, a related tool to build logic models from discovery mass-spectrometry based Phosphoproteomic data. PHONEMeS is described in this paper:  |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **CDA Terfve, E Wilkes, P Casado, P R Cutillas, J Saez-Rodriguez.** [Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data](http://www.nature.com/ncomms/2015/150910/ncomms9033/full/ncomms9033.html). _Nature Communications, 2015, 6:8033_ [PDF](http://www.nature.com/ncomms/2015/150910/ncomms9033/pdf/ncomms9033.pdf). You can also visit [PHONEMeS dedicated webpage](https://saezlab.github.io/PHONEMeS/) |


## CellNOpt Implementations

### _CellNOptR (R packages)_
A series of packages are available in R. The core CellNOpt is available on BioConductor web site: CellNOptR, revision 1.4.0. Newest and oldest version are also available in our [Downloads](http://www.ebi.ac.uk/saezrodriguez/cno/downloads.html) page.

[CellNOptR](http://www.bioconductor.org/packages/release/bioc/html/CellNOptR.html) contains the core functions as well as the boolean and steady states version. It implements the workflow described in Saez-Rodriguez et al Mol Sys Bio 2009, with extended capabilities for multiple time points.

[CNORdt](http://www.bioconductor.org/packages/release/bioc/html/CNORdt.html) is an extension that allows to train a Boolean model agains time-courses of data.

[CNORfuzzy](http://www.bioconductor.org/packages/release/bioc/html/CNORfuzzy.html) is an extension to CellNOptR that allows to handle continous values, using constrained fuzzy logic, as described in Morris et al Plos Comp Bio 2011.

[CNORode](http://www.bioconductor.org/packages/release/bioc/html/CNORode.html) is an ODE add-on to CellNOptR. It is based on the method of (Wittmann et al BMC Sys Bio 2009), also implemented in the tool Odefy (Krusiek et al BMC Bioinf 2010).

[CNORfeeder](http://www.bioconductor.org/packages/release/bioc/html/CNORfeeder.html) is an add-on to CellNOptR that permits to extend a network derived from literature with links derived in a strictly data-driven way and supported by protein-protein interactions as described in (Eduati et al Bioinformatics 2012).

<img src="/cellnopt/public/indexImpl.png" alt="Example Implementations">


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


## Documentation

### _Manual and Tutorial of the R packages_
The R packages are self documented. Tutorials and manual are provided on the bioconductor site of each package. Here below are direct links to the Bioconductor vignettes:

* [CellNOptR vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/CellNOptR/inst/doc/CellNOptR-vignette.pdf)
* [CNORdt vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/CNORdt/inst/doc/CNORdt-vignette.pdf)
* [CNORfuzzy vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/CNORfuzzy/inst/doc/CNORfuzzy-vignette.pdf)
* [CNORode vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/CNORode/inst/doc/CNORode-vignette.pdf)
* [CNORfeeder vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/CNORfeeder/inst/doc/CNORfeeder-vignette.pdf)

Some extra materials and courses about the formats used can be found in the [CNODocs](). Besides, the following link provides a tutorial given at [In Silico Systems Biology, 2013](http://nbviewer.jupyter.org/github/saezlab/cellnopt/blob/gh-pages/public/tutorial_wtac_2013.pdf). The following link provides also a [CytoCopteR tutorial](http://nbviewer.jupyter.org/github/saezlab/cellnopt/blob/gh-pages/public/CytocopterManual.pdf).


## References
**Main reference describing CellNOpt, which can be used to cite it:** 

+ C Terfve, T Cokelaer, A MacNamara, D Henriques, E Goncalves, MK Morris, M van Iersel, DA Lauffenburger, J Saez-Rodriguez. [CellNOptR: a flexible toolkit to train protein signaling networks to data using multiple logic formalisms](http://www.biomedcentral.com/1752-0509/6/133/abstract). _BMC Systems Biology, 2012, 6:133_ [PDF](http://www.biomedcentral.com/content/pdf/1752-0509-6-133.pdf)


**An overview of the different model formalisms available in CellNOpt:** 

+ A. MacNamara, C. Terfve, D. Henriques, B. Peñalver Bernabé, J. Saez-Rodriguez. [State-time spectrum of signal transduction logic models](http://iopscience.iop.org/1478-3975/9/4/045003). _Physical Biology, 9 045003, 2012._ [Reprint](http://iopscience.iop.org/1478-3975/9/4/045003/pdf/1478-3975_9_4_045003.pdf).


**For the core workflow implemented for Boolean modeling:** 

+ J. Saez-Rodriguez*, L. G. Alexopoulos*, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. [Discrete logic modelling as a means to link protein signalling networks with functional analysis of mammalian signal transduction](http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2824489&tool=pmcentrez&rendertype=abstract). Molecular Systems Biology, 5:331, 2009._[*these authors contributed equally to this work]._


**For the constrained Fuzzy Logic implementation:** 

+ M. K. Morris, J. Saez-Rodriguez, D. Clarke, P. K. Sorger, D. A. Lauffenburger. [Training Signaling Pathway Maps to Biochemical Data with Constrained Fuzzy Logic: Quantitative Analysis of Liver Cell Responses to Inflammatory Stimuli](http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3048376&tool=pmcentrez&rendertype=abstract). _PLoS Comp. Biol., 7(3): e1001099, 2011._


**For the CNORFeeder package:** 

+ F. Eduati, J. De Las Rivas, B. Di Camillo, G. Toffolo, J. Saez-Rodriguez. [Integrating literature-constrained and data-driven inference of signalling networks](http://bioinformatics.oxfordjournals.org/content/28/18/2311) _Bioinformatics, 2012 , 28 (18)_


**For the MATLAB toolbox:** 

+ M. K. Morris, I. Melas, J. Saez-Rodriguez. [Construction of cell type-specific logic models of signaling networks using CellNetOptimizer](http://www.ebi.ac.uk/saezrodriguez/files/Morrisetal2011.pdf). _to appear in Methods in Molecular Biology:Computational Toxicology, Ed. B. Reisfeld and A. Mayeno, Humana Press._ [Preprint](http://www.ebi.ac.uk/saezrodriguez/files/Morrisetal2011.pdf)

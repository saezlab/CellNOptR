---
layout: default
title: Home

---


# CellNOpt

## Overview
**CellNOpt** (from CellNetOptimizer; a.k.a. CNO) is a software used for creating logic-based models of signal transduction networks using different logic formalisms (Boolean, Fuzzy, or differential equations). CellNOpt uses information on signaling pathways encoded as a Prior Knowledge Network, and trains it against high-throughput biochemical data to create cell-specific models.

CellNOpt is freely available under GPL license in R and Matlab languages. It can be also accessed through a python wrapper, and a Cytoscape plugin called [CytoCopter](http://www.cellnopt.org/7_CytoCopter/) provides a graphical user interface.

<img src="{{ site.url }}{{ site.baseurl }}public/index1.png" alt="Example result">




## References
Please use this reference to cite CellNOpt:

> E Gjerga, P Trairatphisan, A Gabor, H Koch, C Chevalier, F Ceccarelli, A Dugourd, A Mitsos, J Saez-Rodriguez, [Converting networks to predictive logic models from perturbation signalling data with CellNOpt](https://academic.oup.com/bioinformatics/article/36/16/4523/5855133). _Bioinformatics_, Volume 36, Issue 16, 15 August 2020, Pages 4523–4524,[PDF](https://academic.oup.com/bioinformatics/article-pdf/36/16/4523/33965427/btaa561.pdf), (open access version on [BioRXiv](https://www.biorxiv.org/content/10.1101/2020.03.04.976852v1))

> C Terfve, T Cokelaer, A MacNamara, D Henriques, E Goncalves, MK Morris, M van Iersel, DA Lauffenburger, J Saez-Rodriguez. [CellNOptR: a flexible toolkit to train protein signaling networks to data using multiple logic formalisms](http://www.biomedcentral.com/1752-0509/6/133/abstract). _BMC Systems Biology_, 2012, **6**:133 [PDF](http://www.biomedcentral.com/content/pdf/1752-0509-6-133.pdf)



```

@article{Gjerga2020cellnoptr,
    author = {Gjerga, Enio and Trairatphisan, Panuwat and Gabor, Attila and Koch, Hermann and Chevalier, Celine and Ceccarelli, Franceco and Dugourd, Aurelien and Mitsos, Alexander and Saez-Rodriguez, Julio},
    title = "{Converting networks to predictive logic models from perturbation signalling data with CellNOpt}",
    journal = {Bioinformatics},
    volume = {36},
    number = {16},
    pages = {4523-4524},
    year = {2020},
    month = {06},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btaa561},
    url = {https://doi.org/10.1093/bioinformatics/btaa561},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/36/16/4523/33965427/btaa561.pdf}
}
```

We have also developed PHONEMeS, a related tool to build logic models from discovery mass-spectrometry based Phosphoproteomic data. PHONEMeS is described in this paper:
 > CDA Terfve, E Wilkes, P Casado, P R Cutillas, J Saez-Rodriguez. [Large-scale models of signal propagation in human cells derived from discovery phosphoproteomic data](http://www.nature.com/ncomms/2015/150910/ncomms9033/full/ncomms9033.html). _Nature Communications, 2015, 6:8033_ [PDF](http://www.nature.com/ncomms/2015/150910/ncomms9033/pdf/ncomms9033.pdf). You can also visit [PHONEMeS dedicated webpage](https://saezlab.github.io/PHONEMeS/) |


## CellNOpt Implementations

### _CellNOptR (R packages)_
A series of packages are available in R. The core CellNOpt is available on BioConductor web site: CellNOptR, revision 1.4.0. The most recent updates from [_Gjerga, Trairatphisan, Gabor et al. 2020_](https://academic.oup.com/bioinformatics/article/36/16/4523/5855133) are available [here](https://github.com/saezlab/cellnopt).

[CellNOptR](http://www.bioconductor.org/packages/release/bioc/html/CellNOptR.html) contains the core functions as well as the boolean and steady states version. It implements the workflow described in Saez-Rodriguez et al Mol Sys Bio 2009, with extended capabilities for multiple time points.

[CNORdt](http://www.bioconductor.org/packages/release/bioc/html/CNORdt.html) is an extension that allows to train a Boolean model agains time-courses of data.

[CNORfuzzy](http://www.bioconductor.org/packages/release/bioc/html/CNORfuzzy.html) is an extension to CellNOptR that allows to handle continous values, using constrained fuzzy logic, as described in Morris et al Plos Comp Bio 2011.

[CNORprob](https://github.com/saezlab/CNORprob) is a probabilistic logic variant of CellNOpt which allows for quantitative optimisation of logical network for (quasi-)steady-state data as described in Gjerga, Trairatphisan, Gabor et al 2020.

[CNORode](http://www.bioconductor.org/packages/release/bioc/html/CNORode.html) is an ODE add-on to CellNOptR. It is based on the method of (Wittmann et al BMC Sys Bio 2009), also implemented in the tool Odefy (Krusiek et al BMC Bioinf 2010).

[CNORfeeder](http://www.bioconductor.org/packages/release/bioc/html/CNORfeeder.html) is an add-on to CellNOptR that permits to extend a network derived from literature with links derived in a strictly data-driven way and supported by protein-protein interactions as described in (Eduati et al Bioinformatics 2012). The most recent version of CNORfeeder, which can also be applied to timecourse data with a logic ordinary differential equations (ODE) formalism can be found [here](https://github.com/saezlab/CellNOpt-Feeder) (Gjerga, Trairatphisan, Gabor et al 2020).

[CellNOptR-MaBoSS](https://github.com/saezlab/CellNOptR-MaBoSS) is a method for training boolean logic models of signalling networks using prior knowledge networks and perturbation data with a stochastic simulator (Gjerga, Trairatphisan, Gabor et al 2020).

[ShinyCNOR](https://saezlab.shinyapps.io/shinycnor/) allows using CellNOptR in an interactive way without coding (Gjerga, Trairatphisan, Gabor et al 2020).

<img src="{{ site.url }}{{ site.baseurl }}public/indexImpl.png" alt="Example Implementations">

Additional features to the CellNOpt modelling family have been introduced: **CellNOpt-ILP**, **CellNOpt-MaBoSS**, **CNORprob**, **Dynamic-Feeder**, **Post-hoc analysis** and the **Shiny application** (see figure below). These new features have been described in detail in the _[Gjerga, Trairatphisan, Gabor et al.]((https://www.biorxiv.org/content/10.1101/2020.03.04.976852v1))_ study.

<img src="{{ site.url }}{{ site.baseurl }}public/CNOv2.jpeg" alt="New features of CellNOpt">


### _MATLAB_
Some features of CellNOpt are also available as a MATLAB toolbox, along with the toolbox Q2LM to analyze models, [here](https://github.com/saezlab/MATLAB-CellNOpt)

### _Python_
A Python package called cellnopt.wrapper provides a python interace to the R packages (CellNOptR, CNORode and CNORfuzzy). It uses rpy2 and is available on [Pypi](http://pypi.python.org/pypi/cellnopt.wrapper/). For more details see the sphinx documentation in the ./doc directory after [downloading](https://pypi.python.org/packages/19/3b/d681c432cebbe482c472eb211a6e4de5fc3e444918be4f173335da769762/cellnopt.wrapper-1.0.5.tar.gz#md5=2828b8498acd4a49e7ec2f9fe19aa551) the wrapper. In addition a pure Python version is developed on [github](http://github.com/cellnopt/cellnopt).

### _Cytoscape Plugin (CytoCopter)_
CytoCopteR is a Graphical User Interface designed as a Cytoscape plugin. It provides an interface to CellNOptR using Rserve. More information is available on [CytoCopter](http://www.cellnopt.org/7_CytoCopter/) page.


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

Some extra materials and courses about the formats used can be found in the [CNODocs](http://www.cellnopt.org/6_CNODocs/). Besides, the following link provides a tutorial given at [In Silico Systems Biology, 2013](http://nbviewer.jupyter.org/github/saezlab/cellnopt/blob/gh-pages/public/tutorial_wtac_2013.pdf). The following link provides also a [CytoCopteR tutorial](http://nbviewer.jupyter.org/github/saezlab/cellnopt/blob/gh-pages/public/CytocopterManual.pdf).


## Main Reference
**Main reference describing CellNOpt, which can be used to cite it:** 

+ C Terfve, T Cokelaer, A MacNamara, D Henriques, E Goncalves, MK Morris, M van Iersel, DA Lauffenburger, J Saez-Rodriguez. [CellNOptR: a flexible toolkit to train protein signaling networks to data using multiple logic formalisms](http://www.biomedcentral.com/1752-0509/6/133/abstract). _BMC Systems Biology, 2012, 6:133_ [PDF](http://www.biomedcentral.com/content/pdf/1752-0509-6-133.pdf)

+ E Gjerga, P Trairatphisan, A Gabor, H Koch, C Chevalier, F Ceccarelli, A Dugourd, A Mitsos, J Saez-Rodriguez, [Converting networks to predictive logic models from perturbation signalling data with CellNOpt](https://academic.oup.com/bioinformatics/article/36/16/4523/5855133). _Bioinformatics_, Volume 36, Issue 16, 15 August 2020, Pages 4523–4524.

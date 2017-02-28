---
layout: page
title: CNODocs
---

**CNODocs contains miscellaneous information that could be useful for users and developers of CellNOpt Software. It is not exhaustive. Other sources of information are contained in the Manuals of each of the packages.**

*Section author: Thomas Cokelaer*

# Tutorials

##  I. Manipulating MIDAS files
First, load the library:

```R
library(CellNOptR)
```

First, let us get some data from ([cellnopt.data](http://saezlab.github.io/cellnopt/5_Models%20and%20Documentation/)), load it and convert it to a CNOlist data structure, which is the common data structure used in CellNOptR.:

```R
cnolist = CNOlist(CNOdata("MD-ToyMMB.csv"))
```

**Note:** The CNOdata function search for the file on the cellnopt.data repository. If not found, it looks on your local directory.

Then, we can look at the content of the cnolist variable using the plotting function:

```R
plot(cnolist)
```

<img src="/cellnopt/public/Tutorials1.png" alt="Example Tutorials 1">

The columns contains each species that has been measured. Each row correspond to different conditions summarized in the last column (cues).

You can get some information about the cnolist variable by printing it:

```R
>>> print(cnolist)
class: CNOlist
cues: EGF TNFa Raf PI3K
inhibitors: Raf PI3K
stimuli: EGF TNFa
timepoints: 0 10
signals: Akt Hsp27 NFkB Erk p90RSK Jnk cJun
variances: Akt Hsp27 NFkB Erk p90RSK Jnk cJun
--
To see the values of any data contained in this instance, just use the
method (e.g., getCues(cnolist), getSignals(cnolist), getVariances(cnolist), ...
```

And access to one of the field (e.g. signals) by typing:

```R
>>> cnolist@signals
```


##  II. Manipulating PKN (prior knowledge network)

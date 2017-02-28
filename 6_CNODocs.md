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
Let us first load a model that should be in SIF format (3 columns CSV like):

```R
library(CellNOptR)
pknmodel = readSIF(CNOdata("PKN-ToyMMB.sif"))
```

and plot it:

```R
plotModel(pknmodel, output="PNG")
```

<img src="/cellnopt/public/Tutorials2.png" alt="Example Tutorials 2">

The red edges indicates an inhibitor whereas other black edges are normal links. As we will see in the next section, you can add a cnolist as an argument to colorize this plot.

Using the previous CNOlist data structure, you can color the nodes according to their functions:

```R
cnolist = CNOlist(CNOdata("MD-ToyMMB.csv"))
plotModel(pknmodel, cnolist)
```

Red correspond to inhibitors, green to ligands and blue to measurements.


##  III. Boolean logic (one steady state)
Let us now optimise a MIDAS data set to a PKN network. First, let us load the data provided within CellNOpt:

```R
library(CellNOptR)              # load the library
data(CNOlistToy, package="CellNOptR")                # load the data
cnolist = CNOlist(CNOlistToy)   # convert to a CNOlist instance
data(ToyModel, package="CellNOptR")                  # load the model
pknmodel = ToyModel             # rename it
```

You can visualize the Network (with data information) as follows:

```R
plotModel(pknmodel, cnolist)
```

<img src="/cellnopt/public/Tutorials3.png" alt="Example Tutorials 3">

Preprocess the data before the optimisation:

```R
model = preprocessing(cnolist, pknmodel, compression=TRUE)
```

**Note:** The preprocessing performs 3 tasks. First, it removes nodes that are not connected to any cues (green nodes) or readouts (blue nodes). Second, it compresses nodes that are not cues, readouts or inhibitors (red nodes) except those that have multiple inputs AND ouputs (e.g. Mek node in this case). Finally, it expands the nodes that have multiple inputs by adding “and” gates.

As an exercice, you can plot the new processed model.

Now, let us optimise the model against the data with a Genetic Algorithm. There are many optimisation method available in the literature but a GA is good enough for this toy problem. In the following statement with use the gaBinaryT1 function that performs a boolean optimisation on the first time point only:

```R
opt = gaBinaryT1(cnolist, model, verbose=TRUE)
```

and look at the scores over time:

```R
plotFit(opt)
```

<img src="/cellnopt/public/Tutorials4.png" alt="Example Tutorials 4">

Since it has converged, we can now plot a figure that decomposes the fit of the best model against the data:

```R
cutAndPlot(cnolist, model, bStrings=list(opt$bString))
```

In the resulting figure we can see that the specy NFkB has a high mismatch, which is indicated by the colored boxes

<img src="/cellnopt/public/Tutorials5.png" alt="Example Tutorials 5">


##  IV. Boolean logic (two steady states)

Finally, we can repeat the previous procedure by including an analysis on a second time point. This can be done with a new set of data that includes 2 time points.

```R
library(CellNOptR)              # load the library
data(CNOlistToy2)
cnolist = CNOlist(CNOlistToy2)
data(ToyModel2, package="CellNOptR")
pknmodel = ToyModel2

# the preprocessing
model = preprocessing(cnolist, pknmodel, compression=TRUE)
optT1 = gaBinaryT1(cnolist, model, verbose=FALSE)
optT2 = gaBinaryTN(cnolist, model, verbose=FALSE, bStrings=list(optT1$bString))
```

As before, we can look at the first time points results using:

```R
# looking at the results
cutAndPlot(cnolist, model, bStrings=list(optT1$bString))
```
<img src="/cellnopt/public/Tutorials6.png" alt="Example Tutorials 6">

The second optimisation results can be shown using the same function but providing the two optimisation bitstrings:

```R
# looking at the results for 2 time points
cutAndPlot(cnolist, model, bStrings=list(optT1$bString, optT2$bString))
```

<img src="/cellnopt/public/Tutorials7.png" alt="Example Tutorials 7">

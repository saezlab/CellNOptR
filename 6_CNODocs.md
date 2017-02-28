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


# Formats

##  I. SIF Format
*Section author: Thomas Cokelaer, 2013*

The SIF format (simple interaction format) is a cytoscape compatible format, which is just a space separated value format. It is convenient for building a graph from a list of interactions. It also makes it easy to combine different interaction sets into a larger network, or add new interactions to an existing data set.

The advantage resides in its simplicity:

```R
specy1 1 specy2
specy1 1 specy3
```

with the disadvantage that this format does not include any layout information. Each line in a SIF file represent an interaction between a source and one or more target nodes:

```R
nodeA relationship nodeB
nodeC relationship nodeA
nodeD relationship nodeE nodeF nodeB
```

The SIF format used in CellNOpt is actually a subset of the official SIF format:

```R
#. Generally there is only one target
#. The relationship can be only the number 1 or -1 (for inhibitor)
#. Duplicate entries are not ignored.
```

Delimiters can be spaces or tabs (mixed).

Duplicated entries are ignored. Multiple edges between the same nodes must have different edge types.

Self loops are also possible:

```R
A 1 A
```

In cytoscape, the relationship type can be any string. However, in CellNOptR we are limited to the values 1 and -1 that correspond to activation or inhibition.


##  II. MIDAS Format
*Section author: Thomas Cokelaer, 2013*

This document describes briefly the MIDAS (Minimum Information for Data Analysis in Systems Biology) format that is used in CellNOpt software.

### 2.2.1. Format
MIDAS files are CSV files (comma separated). The content is defined in the first line of the file that constitutes the header (only 1 line). In the header, colums can take two forms:

```R
XX:Specy,
XX:userword:Specy
```

where XX is a 2-letter word prefix that describes the column content (see table below for valid word) and Specy is the name of the column. The userword is optional (see later).

MIDAS files include the concept of cues, signals, and responses (Gaudet et al., 2005):

* cues are biological perturbations to a system (such as the addition of extracellular ligands)
* signals represent the activities of proteins or other biomolecules involved in transducing biological information (activation of an intracellular kinase, for example),
* responses also called readouts represent phenotypic changes such as proliferation, cell death or cytokine release.

The column headers in a MIDAS files may contain a second (userword) level of identification (e.g. headers for columns describing various cytokine treatments might begin with “TR:Cytokine”). When present, these secondary identifiers allow Software (e.g., DataRail’s importer) to identify automatically the dimensions of a new compendium.

| Code | Description | handled in CellNOptR |
| ---- | ----------- | -------------------- |
| ID | Identifiers |  |
| TR | Treatment | yes |
| DA | Data acquisition | yes |
| DV | Data value | yes |

Example:

| TR:mock:CellLine | TR:EGF | TR:TNFa | TR:PI3Ki | DA:Akt | DA:Hsp27 | DV:Akt | DV:Hsp27 |
| ---------------- | ------ | ------- | -------- | ------ | -------- | ------ | -------- |
| 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| 1 | 1 | 0 | 0 | 10 | 10 | 1 | 0.2 |
| 1 | 0 | 1 | 0 | 10 | 10 | 1 | 0.5 |

Each value is separated by a comma and you could have space, tabs between commas. So, the final format could be as follows:

```R
TR:mock:CellLine, TR:EGF, TR:TNFa, TR:PI3Ki, DA:Akt, DA:Hsp27, DV:Akt, DV:Hsp27
1,1,0,0, 0,0,   0,0
1,0,1,0, 0,0,   0,0
1,1,0,0, 10,10, 0.82,0.7
1,0,1,0, 10,10, 0.91,0.7
```

Let us explain the header:
* The first row is the header describing the content of each colum.
* commas separate all fields.
* Each fields must starts with one of the valid code followed by a column (e.g., TR: or DA:).
* Subfields are possible: TR:whatever:EGF
* Special fields such as CellLine, NOCYTO and NOINHIB are ignored.
* The number of DA and DV must be equal except if you use the special name DA:ALL (see later).
* Inhibitors are coded by adding the letter i after the name (e.g., TR:PI3Ki)
**Warning:** do we have special cases of name ending with the letter i ?

The data above is made of rows that length is as long as the header. Fields may be empty, which is not the case here. If so, software should replace the value by (e.g., NA in R language) and cope with it.

Each row represents a given treatment at a given time. Time are coded with the DA code. Values are coded within the **DV** columns. Let us look at the 2 first rows. The time is 0. The next two other rows are coded for the time 10. The treatements (3 first colums) are found at the different time.

In MIDAS file, data should be ordered by time although some software may deal with it.

### 2.2.2. Filename issue
From the *Reference* below:

```R
MIDAS file has a unique identifier (UID) composed of the following fields:
(i) a two-letter data/file-type code (e.g., PDfor Primary Data, MD for
multiplex data), (ii) a three-letter creator code (typically initials),
(iii) an identification number of arbitrary length that is unique across
the entire system, and (iv) a free-text suffix that serves as a mnemonic
to improve human readability. For example, the primary data discussed in
the text might be tagged MD-LGA-11111-CytoInh17phFI-BLK
```

In practice, only a few files are coded that way. One reason is that the UID tag is hardly used. Another inconsistency is that dashes are not used or replaced by ____. Besides, many files contain the word Data. Finally, the name tag (e.g. LGA above) is not good practice because public file should give the feeling they belong to everybody. However, one consistency is the extension being **.csv**.

*Reference*: **J. Saez-Rodriguez, A. Goldsipe, J. Muhlich, L. Alexopoulos, B. Millard, D. A. Lauffenburger, P. K. Sorger**, *Flexible Informatics for Linking Experimental Data to Mathematical Models via DataRail.* Bioinformatics, 24:6, 840-847 (2008). [Citations](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btn018)

### 2.2.3. Proposal for filename convention

* Do not use the DATA/Data word. Instead start all files with MD- and use the extension .csv
* Separate the names that describe your data with dashes.
* Underscore could be use internally to refine a name
* MD must be capitalised, other names can use any convention but we recomment polish convention (e.g., capitalize words)

MD-Tag1-Tag2.csv

MD indicates that this is a MIDAS file so no need to set Data in the filename anymore. Tag1 is a general description tag (containing _ possibly) and Tag2 is a variant of Tag1. For instance, Tag1 could be Toy and Tag2 a name to differentiate different Toy data sets.

Correct:

```R
MD-Toy.csv
MD-Toy-variant1.csv
MD-LiverDream.csv
MD-LiverDREAM.csv
```

[DataRail](https://sites.google.com/site/saezrodriguez/software/datarail)

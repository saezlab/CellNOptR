---
layout: page
title: Usage
---

Running PHONEMeS requires the *[igraph package](http://igraph.org/r/)* to be installed, while in addition for running the [CNORode package](https://bioconductor.org/packages/release/bioc/html/CNORode.html) the *[MEIGO optimization toolbox](https://www.bioconductor.org/packages/release/bioc/html/MEIGOR.html)* needs also to be installed first.

# CellNOpt workflow

## I. Data Input

CellNOpt method derives a logic-based model from a *prior knowledge network* (PKN) and trains it against perturbation measurements. The PKN is typically inferred from information from literature or expert knowledge and then tests this prior knowledge if it is compatible with the measurements in a specific perturbation context. All CNO formalisms, typically take two kind of inputs: **i)**A SIF (Simple Interaction File) text format which containing the *PKN* with the signed directed interactions happening between sites and which can be manipulated in *Cytoscape* and **ii)**A MIDAS (Minimum Information for Data Analysis in System Biology) containing data which are representing measurements of sites upon various stimuli and inhibitory perturbations.

## II. Network Processing

To prepare the prior knowledge network for training, we first apply two pre-processing steps: **i)**In *compression* we reduce and simplify the network by removing species that are neither measured nor perturbed without impairing tha logical consistency of the network. **ii)**In *expansion* all the logical operations from each interactions are added.

## III. Training

Training involves the identification of the sub-models of the processed $PKN$ which better explains and fits the data by minimizing the mean squared deviation between model prediction and the actual measurements. We also add a size penalty to the objective function as the sum of the number of inputs of each edge of the prior knowledge scaffold network, normalized by the total number of inputs across all edges. $ \theta $ Test test test $x$, $y_1$, $ \Theta $

## IV. Report

Networks inferred are written to a file and plotted, while information about the optimization during the training step are summarized in a central HTML report hyperlinked to the different plots.

# Examples

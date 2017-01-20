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

Training involves the identification of the sub-models of the processed $PKN$ which better explains and fits the data by minimizing the mean squared deviation $ \theta_f $ between model prediction $ B^{M} $ and the actual measurements $ B^{E} $ for the $ m $ readouts, $ n $ time-points and $ s $ experimental conditions weighted by the total number of data points $ n_{g} $. 

$$ \theta_f(P) = \frac{1}{n_{g}} \sum_{k=1}^{s} \sum_{l=1}^{m} \sum_{t=1}^{n} (B_{k,l,t}^{M}-B_{k,l,t}^{E})^{2} $$

$ P $ on this case, represent one of the solutions or the string of bits indicating whether an edge of the compressed network is included in the model or not.

We also add a size penalty $ \theta_s $ to the objective function as the sum of the number of inputs $ \nu_{e} $ of each edge $ e $ in $ P $ normalized by the total number of inputs across all edges ($$ \nu _{e}^{s} = \sum_{e=1}^{r} $$) where $ r $ represents the length of $ P $.

$$ \theta_s(P) = \frac{1}{\nu _{e}^{s}} \sum_{e=1}^{r}\nu _{e}P_{e} $$

## IV. Report

Networks inferred are written to a file and plotted, while information about the optimization during the training step are summarized in a central HTML report hyperlinked to the different plots.

# Examples

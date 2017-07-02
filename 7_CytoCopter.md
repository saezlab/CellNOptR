---
layout: page
title: CytoCopter
---

CytoCopteR provides an intuitive and easy to learn graphical user interface (GUI) to [CellNOptR](https://saezlab.github.io/CellNOptR/7_CytoCopter/)

This results in a point and click interface where users can run the same steps as they would using an R script without having to actually write any code. Given that this is a front-end to the [R](https://www.r-project.org/) algorithms, consistency is ensured between the results obtained through the GUI and those obtained through the corresponding scripts.

# Download

CytoCopteR can be downloaded from the cytoscape App store ([Download link here](http://apps.cytoscape.org/apps/cytocopter)). It can also be directly installed from cytoscape  by following these simple steps:

Go to **Apps** and then **App Manager...** and finally search for CytoCopteR in the search bar and press **Install**.

<img src="/CellNOptR/public/cytocopter_1.png" alt="Installing Cytocopter">


# Usage

## Load Model file

After CytoCopteR is installed, we can now import the data files through the import function in Cytoscape. We first import the network [ToyModelPB.sif](http://nbviewer.jupyter.org/github/saezlab/CellNOptR/blob/gh-pages/public/ToyModelPB.sif) through the import function in Cytoscape.

## Preprocess experimental data and network

At this point we already have the network imported, now we can use our experimental data and see how it maps into the network. In other words visualise which nodes are measured, inhibited and stimulated. For more details, please consult [*Terfve et al.*](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-6-133).

  * Select *CytoCopter* panel on the left side panel
  
  * On the network dropdown box, we have our imported network
  
  * Click on the data text fiels and a window will pop-up to browse the previosly downloaded experimental ([MIDAS file](http://nbviewer.jupyter.org/github/saezlab/CellNOptR/blob/gh-pages/public/ToyModelPB.csv))
  
  * After the data is selected, press the *Preprocess* button. **This procedure may take a few minutes; in particular if it is the first time all the necessary R packages will be installed**
  
  * The *Preprocess* function will annotate automatically the network: Green nodes are stimulated species; Red nodes are inhibited species; Blue nodes are measured species; Grey nodes with dashed borders can be removed to simplify the network; Blue nodes with red borders are measured and inhibited nodes; White nodes are not measured or perturbed that cannot be simplified.
  
<img src="/CellNOptR/public/cytocopter_4.png" alt="Preprocessing">

## Optimising network with experimental data

In this step, we will use the optimise features of CellNOptR R package to optimise the topology of the network against previously loaded experimental data. For more details, please consult [*Terfve et al.*](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-6-133).

 * Select a formalism to model the network. Select the Boolean Formalism. Boolean formalism will consider that the states of each node will be either active or inactive. **Because the Boolean formalismrequires anly two time points and the given data-set has 16, the first time point is automatically set to 0 and the second one is selected by the user**
 
 * Select time point 10 from the *Time point* dropdown box
 
 * After selecting the time point, press the optimise button. **Note this will take several minutes. This optimisation is done using a specific genetic algorithm and its parameters can be configured in the Configuration panels, such as maxmum time allowed for the algorithm to run in the Max time parameter**.
 
 * The optimisation procedure outputs an optimised network and a plot displaying how well the simlated models fit with the given experimental data. The colour gradient quantifies the error levels.
 
<img src="/CellNOptR/public/cytocopter_5.png" alt="Optimisation">

 


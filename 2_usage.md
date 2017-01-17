---
layout: page
title: Usage
---

CellNOpt method derives a logic-based model from a *prior knowledge network* (PKN) and trains it against perturbation measurements. The PKN is typically inferred from information from literature or expert knowledge and then tests this prior knowledge if it is compatible with the measurements in a specific perturbation context. All CNO formalisms, typically take two kind of inputs: **i)**A SIF (Simple Interaction File) text format which containing the *PKN* with the signed directed interactions happening between sites and which can be manipulated in *Cytoscape* and **ii)**A MIDAS (Minimum Information for Data Analysis in System Biology) containing data which are representing measurements of sites upon various stimuli and inhibitory perturbations.

# Setting up for SteadyCom

##### Authors: Arya Narnapatti, Namasi G
##### In case of queries, suggestions, further contributions, contact at: arya.narnapatti@students.iiserpune.ac.in
<br>

## Introduction
SteadyCom is an algorithm used to analyze the metabolic behaviors of co-cultures and microbiomes in a community steady state. It assumes the existence of a time-averaged, population steady-state for a stable microbial community which in turn implies a time-averaged, constant growth rate across all members. It also takes into account the biomass abundances of members of the community as opposed to being a simple extension of Flux Balance Analysis (FBA) to communities, by using a multicompartment model. More about the SteadyCom algorithm and its derivation can be found at [ref[1]](#references). SteadyCom.m is a function integrated into the COBRA toolbox [<sup>2</sup>](#references) that allows one to analyze a multi-organism COBRA model using the SteadyCom algorithm.


This document serves to familiarise users with the requirements of a multi-organism COBRA model for SteadyCom and how to build one. It also serves to act as a repository for methods to troubleshoot issues faced while running SteadyCom on COBRA. For more on using SteadyCom to analyze microbial communities, view [ref[3]](#references). 

We also hope to collate more information to add to this document and would encourage other users to add to it as well.
<br>

## The Multi-species Model
*This section talks more about the structure of the joint model required by SteadyCom, and the motivation behind this document.*
 
As an input, SteadyCom takes in a community COBRA model structure that can be created using the *‘createMultipleSpeciesModel’</i>* function on COBRA. This model incorporates a common community compartment (also called lumen) that allows for the exchange of metabolites between the different individual species. Instead of merging the extracellular compartments of all the species in the consortium, this community compartment is made in addition to having separate extracellular compartments for each species. This makes it easier to monitor the net flux of metabolites through the membranes of each species in the joint model using **“shuttle reactions”** between the community compartment and the individual extracellular compartments[<sup>4</sup>](#references). 

The joint model contains a stoichiometric matrix which contains all the reactions of each individual model along with the **shuttle reactions** (between the community compartment and the extracellular compartment of each organism) and the **community exchange reactions** (capturing the net flux of all extracellular metabolites into and out of the consortium). In fact, this stoichiometric model will contain blocks that are nothing but the individual stoichiometric matrices of each species, except that the original exchange reactions will now drain out into the lumen and that there will be separate columns for the community exchange reactions (see [ref[4]](#references) to learn more). Apart from this, the joint model will contain other fields such as upper and lower bounds of all reactions, the objective coefficients(c), etc -  similar to an individual genome-scale model (GSM).

SteadyCom requires at least one of two additional structures to be present in the joint model (we will call these **indCom** and **infoCom** in this document). These contain details about the extracellular metabolites, the shuttle reactions, the community exchange reactions, the biomass reactions of each species, and more. These fields are not present in the model output from the *‘createMultipleSpeciesModel’* function but can be obtained from the *‘getMultiSpeciesModel’</i>* function on COBRA.

<p align="center">
  *
 </p>

The inputs required for the *‘createMultipleSpeciesModel’* function that outputs the joint model are the individual Genome-Scale Models of each species in the consortium. To create the joint model things like extracellular metabolites, exchange reactions, compartments in the model, and more need to be recognized. This means that the individual models need a level of standardization to be compatible with the function. For example, at the most basic level the extracellular metabolites of both models must be named similarly. At times a particular model may name sucrose as 
<span style="color:blue">
  'sucr_e'
</span> 
whereas another might name it cpd00076[e]. In such a case, the function will not be able to recognize that these two are the same and hence will classify them as different in the community compartment, which means interactions molded by sucrose will be missed out. 


## References

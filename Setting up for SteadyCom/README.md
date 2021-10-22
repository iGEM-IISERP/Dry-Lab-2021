
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

The inputs required for the *‘createMultipleSpeciesModel’* function that outputs the joint model are the individual Genome-Scale Models of each species in the consortium. To create the joint model things like extracellular metabolites, exchange reactions, compartments in the model, and more need to be recognized. This means that the individual models need a level of standardization to be compatible with the function. For example, at the most basic level the extracellular metabolites of both models must be named similarly. At times a particular model may name sucrose as `'sucr_e'` whereas another might name it `'cpd00076[e]'`. In such a case, the function will not be able to recognize that these two are the same and hence will classify them as different in the community compartment, which means interactions molded by sucrose will be missed out. 

Similarly, exchange reactions, compartments, etc are named or structured differently in different models and can cause issues. We have thus prepared a short checklist that can help users avoid or troubleshoot some of these issues. Starting from polishing the individual models, to setting up the joint model and certain methods that users can try to troubleshoot issues with running SteadyCom. 

*Notes:*<br>
*To begin working with the COBRA toolbox on MATLAB, one has to initialize the toolbox.*<br>
*Models can be loaded using the ‘readCbModel’ function. (Usage at [ref[5]](#references))*<br>
*For example in this document, individual models will be loaded as `model`,  `model1`, `model2`, etc (as the variable name) and the joint model will usually be loaded as `EcCom`.*<br>

## The Checklist

### Preparing the Individual Models 
1.
An easy way to get a quick overview of the individual models or to extract simple information from them without worrying much about the various COBRA functions is to save the model as a .xls/spreadsheet file. This can be done using the *‘writeCbModel’* function: 
```matlab
writeCbModel(model1,'.xls')
```
This will save `model1` as an .xls file on the directory open on MATLAB. We recommend that the user save the individual models as .xls files as they would be helpful in some of the subsequent steps. (Though it isn’t required as the same steps can be done using functions on MATLAB). The .xls file created will have two sheets - one listing all the reactions of the model along with their descriptions, formulas, bounds, etc and the other listing the metabolites of the model along with their descriptions, formulas, compartment, etc. (Note: the same metabolite in the different compartments are listed as different metabolites.)

2.
One should make sure to check the number of extracellular metabolites and exchange reactions in the model. Ideally, there should be a one-to-one mapping where each extracellular metabolite has a corresponding exchange reaction, but in many models, there are extracellular metabolites that do not have the corresponding exchange reactions. (These probably do not cause issues while running FBA optimizing for growth as these metabolites do not need to be produced or taken in to grow but can cause issues while setting bounds for the shuttle reactions in the joint model later). The *‘createMultipleSpeciesModel’* function recognizes exchange reactions as those that contain ‘EX’ in the name and are for extracellular metabolites. In some models, other sink or demand reactions that are not exchange reactions (those which are for metabolites in the inner compartments of the cell) also have ‘EX’ in their names. This can also cause problems as those reactions will be missing from the joint model.

The first thing to make sure is that the number of true exchange reactions equals the number of extracellular metabolites. This can be done using the .xls file. Exchange reactions can be found by searching for ‘EX_’ through the first column of the reaction list. Make sure to keep an eye out for sink or demand reactions containing ‘EX_’, for the aforementioned reasons. It may also be possible that some models contain exchange reactions without ‘EX_’ . To identify such reactions, one can use the *‘findExcRxns’* function that identifies all reactions with just one metabolite in the formula (so this will include all exchange, sink, demand reactions etc.). The reactions here that correspond to extracellular metabolites are considered to be **true exchange reactions**. 

The extracellular metabolites also usually have ‘[e]’ at the end of their names. In some models, they may have ‘_e’ or ‘[extracellular]’, etc. Searching for this string in the first column of the metabolite list can allow one to count the number of extracellular metabolites in the model. Note that this isn’t a perfect method. For example, there may be other metabolites with _e in the name, or certain metabolites might be repeated or the naming system might not be standardized.

When exchange reactions for certain metabolites are missing, they can be manually added to the models. This can be done using the ‘addReaction’ function mentioned in point 5 below.

3.
**Compartments** - To ensure compatibility of the models with reconstructions from the BIGG database[<sup>6</sup>](#references), compartments should be identified correctly, as according to - http://bigg.ucsd.edu/compartments. Further, metabolite names should end with [e], [c], etc. to identify them with their respective compartments and not _e, _c, or other naming conventions. The *‘createMultiSpeciesModel’* function itself has a section to edit the metabolite names to ensure compatibility, but for some models with uncommon conventions, there may be issues. A possible troubleshoot is to rename metabolites and compartments according to the BIGG convention before using the *‘createMultiSpeciesModel’* function. Here is an example where metabolites were identified with ‘_e’ instead of [e],  ‘_c’ instead ‘[c]’, etc .
```matlab
model.mets = regexprep(model.mets, '_([^_])$', '\[$1\]');
```

4.
One of the most common issues faced when creating a multi-species model is the fact that the naming conventions for metabolites in different models are different. Some models have metabolites named according to the modelSEED convention (e.g `’cpd00076’`, `’cpd00214’`, etc) whereas others have metabolites named as`’ sucr_e’`, `’glc__D_e’`, etc. For internal metabolites this isn't an issue, but as mentioned before common extracellular metabolites have to be named in the same way. 
The most reliable way to do this is to select the convention of one of the individual models and then manually go through the extracellular metabolites of the rest of the models and rename them according to the chosen convention. One method of renaming a metabolite name is: 
```matlab
model.mets(findMetIDs(model, ‘oldmetname’)) = {‘newmetname’}
```

5.
In many cases, users may need to modify reconstructed individual species models before integrating them into the joint model. This may be due to missing reactions in the model, strain differences, or gene modifications made for any purpose. Users may wish to add or delete reactions, metabolites, or genes of a specific function to the present model. Those familiar with these basic COBRA Toolbox functions may skip through this section. 
> a. **Adding or deleting reactions:** To add a reaction to a model, the *‘addReaction’* function can be used in two ways: One of the approaches is to type out the entire reaction as a string. One must be careful with using spaces wherever necessary and try to preserve the same string format when using this method. It is, in general, important to preserve the naming conventions of reactions (for example, using ‘EX_’ prefix for exchange reactions) when adding reactions. In the example given below, we see that the function takes the model, reaction name, and a name-value pair of the reaction formula as the input:
> ```matlab
>model = addReaction(model, 'GAPDH', 'reactionFormula', 'g3p[c] + nad[c] + 2 pi[c] -> nadh[c] + h[c] + 13bpg[c]');
>```
>Another approach is to add reactions as a list. Here, the function takes care of string formatting. The inputs are model, reaction name, name-value pair of metabolite list, name-value pair of the stoichiometric coefficients of each of the metabolites respectively, and the reversibility of the reaction.
> ```matlab
>model = addReaction(model, 'GAPDH2','metaboliteList', {'g3p[c]', 'nad[c]', 'pi[c]', '13bpg[c]', 'nadh[c]', 'h[c]' },...
>'stoichCoeffList', [-1; -1; -2; 1; 1; 1], 'reversible', false);
>```
>Note: The addReaction function automatically adds any new metabolites present in the reaction to the list of metabolites. This is why it’s important to make sure spelling and string format are maintained.If any new metabolites are being added, make sure there is a source that they are produced from or an uptake reaction is provided into the extracellular compartment along with necessary transport reactions. For any metabolite that needs to be secreted, make sure there are the necessary transport reactions into the extracellular compartment and an exchange reaction from this compartment that allows the metabolite to be drained out.
>
>To delete a reaction, use the *'removeRxns'* function. A list of reactions can be removed at once by providing all the reaction names inside {}.
> ```matlab
>model = removeRxns(model, {'EX_glc_D[c]', 'EX_glc_D[e]', 'sink_13bpg[c]', 'sink_nad[c]', 'DM_dhap[c]', 'DM_g3p[c]'});
>```
>Note that it is recomended that knockouts be simulated not by removing reactions but by just setting both the upper and lower bounds of the reactions to 0. Reactions should be >removed only when required.

> b. **Adding or deleting metabolites:** As mentioned previously, metabolites are added automatically when reactions that involve these new metabolites are added.
To remove metabolites from the model, use the *‘removeMetabolites()’* function. The last input to this function is true or false which defines whether the user would like to remove all the reactions pertaining to the metabolite or simply remove the metabolite from the reactions wherever present. ‘True’ dictates the reactions to be deleted, while ‘false’ dictates only the metabolite to be removed while preserving the rest of the reaction.
> ```matlab
>model = removeMetabolites(model, {'3pg[c]', '2pg[c]'}, false);
>```

> c. **Adding or deleting genes:** For SteadyCom genes of a model are not a required structure that needs to be present. However, many models have a field with the names of the genes and the reactions associated with them to group together reactions and understand the model easier. 
Genes are added by associating reactions to a gene name. This can be done in the*’ addReaction function’* itself uses a name-value pair `‘geneRule’`. `‘or’` and `‘and’` are used to associate multiple genes to a reaction when the gene rules are complicated.
> ```matlab
>model = addReaction(model, 'PGK', 'geneRule', 'G2 or G3');
>```
>To associate a reaction already present in the model with a gene, or to change the gene associated with a reaction, *’changeGeneAssociation’* is used.
> ```matlab
>model = changeGeneAssociation(model, 'GAPDH', 'G1 and G2'); 
>```
>To delete a gene, it is recommended to simply set all the reactions associated with that gene to flux 0, by setting lowerbound (lb) and upperbound (ub) of all these reactions to 0. This is useful as it can be referred to later and the genes and reactions are not completely lost from the model. 

6.
It is recommended to run FBA on individual models once all the modifications to the models have been completed. This acts as an important verification to make sure no errors are present in the model. A user could try to set different objective functions, such as the growth of biomass, or any other exchange reaction, and compare the fluxes to make sure they provide realistic values. One could try to verify these values by comparing them with experiments done with these strains in a lab (or from literature) to make sure the flux values from the GSM are comparable to the real biomass growth rate of the strain or rate of secretion of certain metabolites.

### Setting up the Co-culture Model 
1.
**Basic Steps to setting up the co-culture model**
> I. **Setting up the nametags** <br>
> This sets a nametag for each individual model that will be used as a prefix for the metabolites and reactions of that model in the joint model.
> ```matlab
>nameTagsModel = {'model1'; 'model2'};
>```
>Instead of model1, model2, one can use nametags like ‘ecoli’, ‘yeast’, etc as well.

> II. ***‘createMultipleSpeciesModel’***
> ```matlab
>EcCom = createMultipleSpeciesModel({model1; model2}, nameTagsModel);
>```
>*(Here model1 and model2 are the names of the variables the individual models are stored as)*

> III. ***‘getMultiSpeciesModelId’***
>*‘getMultiSpeciesModelId’* is used to form the `indCom` and `indoCom` structures inside the joint model
> ```matlab
>[EcCom.infoCom, EcCom.indCom] = getMultiSpeciesModelId(EcCom, nameTagsModel);
>```

> IV. **Setting the individual biomass reactions**
> ```matlab
>rxnBiomass = {'model1Biomassrxnname', 'model2Biomassrxnname}.'; % Identifying the biomass reaction names from each model
>rxnBiomassId = findRxnIDs(EcCom, rxnBiomass); % getting the corresponding rxnID’s 
>EcCom.infoCom.spBm = rxnBiomass; % Storing the biomass reaction names in infoCom
>EcCom.indCom.spBm = rxnBiomassId; % Storing the biomass reaction ID’s in indCom
>```
>Note here that `‘modelxBiomassrxnname'` should be the name of the biomass reactions of `modelx`

2.
The joint model (named `EcCom` in this document) made at the end of these 4 steps will have many of the fields that individual genome-scale models have such as S (the stoichiometric matrix), c (the objective coefficients), csense, lb and ub(the bounds) and more. Two additional structures that are special to the joint model are `indCom` and `infoCom`. (Note that at least one of these two is required for SteadyCom, but the *‘getMultiSpeciesModelId’* function is capable of forming both). Both these structures have subfields that hold information about all the extracellular metabolites (Mcom), shuttle reactions (EXsp), community exchange reactions (EXcom), biomass reactions (spBm), etc. 

infoCom stores information using strings (metabolite names, reaction names, etc) while indCom stores numerical values (reaction IDs, metabolite IDs, etc). 

EXsp holds information about the shuttle reactions of the model. Each column corresponds to an individual species in the joint model and each row to a metabolite present in the extracellular space. The (i, j)<sup>th</sup> entry of EXsp will hold information about the shuttle reactions for the ith metabolite in the j<sup>th</sup> species. The name of the i<sup>th</sup> metabolite can be found from `infoCom.Mcom`, which contains the list of all extracellular metabolites. A zero or blank at the(i, j)<sup>th</sup> position indicates that the i<sup>th</sup> metabolite is not taken up or given out by the j<sup>th</sup> species (basically, the metabolite is absent in the GSM of that species) and hence will not play a role in any interactions of that species.

The number of external metabolites in the individual models should match the number of shuttle reactions for each species in the joint model. The shuttle reactions are named in the format `‘model1IEX_metname[u]tr’` where model1 here is the individual model. 
*Note: The shuttle reactions are named according to the external metabolite names and not the exchange reaction names of the individual models.* 

3.
**Checking the csense** <br>
csense is another structure in the joint model. It is a vector that contains the constraint sense for each row in the stoichiometric matrix. ‘E’ denotes equality, ‘G’ for greater than, and ‘L’ for less than. At times, discrepancies in the csense might cause problems during the SteadyCom run. As default, when all constraint senses are equality, csense should be a character array of length equal to the number of metabolites in the model; and each entry should be E.

4.
**Setting bounds on the community exchange and shuttle reactions to simulate biologically realistic constraints** <br>
While running FBA on individual models, not all exchange reactions are reversible (their lower bound is set to 0). This means that certain extracellular metabolites cannot be taken in by the species, but only given out. The reason for adding these constraints is to make the model more realistic.

However, when the joint model is created through the *’createMultiSpeciesModel’* function, all the shuttle and community exchange reactions are set to be reversible. Hence before running SteadyCom, one needs to constrain these reactions to get a more realistic solution. Ideally, these bounds should be set manually depending on the conditions of growth (media, diet, aerobic/anaerobic). But if the user is not sure about these parameters, a preliminary method would be to set the bounds of the shuttle reactions to be the same as the original individual model used.

 ```matlab
var = EcCom.rxns(nonzeros(EcCom.indCom.EXsp(:,1)));
var = erase(var, 'model1I');
var= erase(var, '[u]tr');
var = append(var, '_e');
var2 = model1.lb(findRxnIDs(model1,var));
EcCom.lb(nonzeros(EcCom.indCom.EXsp(:,1))) = var2;
```

This is an example where shuttle reactions in the joint model `EcCom` and exchange reactions in the individual model `model1` are named according to the extracellular metabolite. 


## Running SteadyCom and Troubleshooting
The basic usage of the SteadyCom function can be accessed by running:
 ```matlab
[sol, result] = SteadyCom(EcCom) 
```

There are additional optional inputs that can be accessed through the `options` structure. This allows users to impose certain additional constraints, choose the specific algorithm they wish to use, and more. When the `options` structure is not input, the default values are taken. This should work for a preliminary run.

**Some possible strategies to troubleshoot when SteadyCom gives an infeasilbe result**<br>
1. <br>
**Changing the algorithm** <br> 
The `options` structure allowed users to select one of three algorithms. By default, this is set to algorithm 1 - a simple guess for bounds followed by Matlab fzero. This is usually the most efficient, however, there are cases when algorithm 1 can face issues, and switching to algorithm 2 or 3 will reach a convergent solution. The algorithm can be changed by setting `options.algorithm` to 2 or 3 

The `options` structure allowed users to select one of three algorithms. By default, this is set to algorithm 1 - a simple guess for bounds followed by Matlab fzero. This is usually the most efficient. However, there are cases when algorithm 1 can face issues, and switching to algorithm 2 or 3 will reach a convergent solution. The algorithm can be changed by setting `options.algorithm` to 2 or 3

For example:
 ```matlab
options.algorithm = 2;
[sol, result] = SteadyCom(EcCom, options)
```

2.

**Using a different solver** <br>
If changing the algorithm doesn't work, changing the solver can also help troubleshoot an infeasible result on SteadyCom. Glpk, gurobi, ibm_cplex are some of the solvers that are compatible with SteadyCom on COBRA.

3.
**Checking for missing sink reactions** <br>
As mentioned before, sometimes sink reactions in the individual models that are in the form ` metname[u]  →  `  (where ‘u’ is an internal compartment of one of the individual GSM’s) do not find their way to the joint model. Sometimes these reactions are required to have flux for the model to grow and hence their absence in the joint model could cause infeasibility when running SteadyCom. These reactions would have to be manually added to the joint model.

4.
If none of these work, one can try building joint models consisting of 2 copies of the same individual species GSM. Building such a model for each species in the consortium and checking if these models are able to show similar results to running FBA on the individual GSM’s could provide insights as to where the issue is arising.  

## Future Plans
As mentioned earlier, we hope to continue to collate more information to add to this document and would encourage other users to add to it as well. Some avenues for future expansion of this document that we believe would be very useful for future users are: <br>

1.
Multi-species models can also integrate interactions with the host (for example in the case of gut microbiota) by including a GSM of the host organism in the community. Incorporating this brings forth its own intricacies and may cause a variety of new problems to troubleshoot. 

2.
The `options`  structure input by SteadyCom. This allows users to impose certain additional constraints, choose the specific algorithm they wish to use, and more. More documentation on this is required.


## References

1. Chan, S. H. J., Simons, M. N., & Maranas, C. D. (2017). SteadyCom: predicting microbial abundances while ensuring community stability. PLoS computational biology, 13(5), e1005539.
2. Laurent Heirendt & Sylvain Arreckx, Thomas Pfau, Sebastian N. Mendoza, Anne Richelle, Almut Heinken, Hulda S. Haraldsdottir, Jacek Wachowiak, Sarah M. Keating, Vanja Vlasov, Stefania Magnusdottir, Chiam Yu Ng, German Preciat, Alise Zagare, Siu H.J. Chan, Maike K. Aurich, Catherine M. Clancy, Jennifer Modamio, John T. Sauls, Alberto Noronha, Aarash Bordbar, Benjamin Cousins, Diana C. El Assal, Luis V. Valcarcel, Inigo Apaolaza, Susan Ghaderi, Masoud Ahookhosh, Marouen Ben Guebila, Andrejs Kostromins, Nicolas Sompairac, Hoai M. Le, Ding Ma, Yuekai Sun, Lin Wang, James T. Yurkovich, Miguel A.P. Oliveira, Phan T. Vuong, Lemmer P. El Assal, Inna Kuperstein, Andrei Zinovyev, H. Scott Hinton, William A. Bryant, Francisco J. Aragon Artacho, Francisco J. Planes, Egils Stalidzans, Alejandro Maass, Santosh Vempala, Michael Hucka, Michael A. Saunders, Costas D. Maranas, Nathan E. Lewis, Thomas Sauter, Bernhard Ø. Palsson, Ines Thiele, Ronan M.T. Fleming, Creation and analysis of biochemical constraint-based models: the COBRA Toolbox v3.0, Nature Protocols, volume 14, pages 639–702, 2019 doi.org/10.1038/s41596-018-0098-2.
3. Chan, S. H. J. Analyze Steady-State Community COBRA Models. (https://opencobra.github.io/cobratoolbox/latest/tutorials/tutorialSteadyCom.html)
4. Klitgord, N., & Segrè, D. (2010). Environments that induce synthetic microbial ecosystems. PLoS computational biology, 6(11), e1001002.
5. https://github.com/opencobra/cobratoolbox/blob/efb2c06be702aaed38349581ba73ed420251992b/src/base/io/readCbModel.m
6. King, Z. A., Lu, J., Dräger, A., Miller, P., Federowicz, S., Lerman, J. A., ... & Lewis, N. E. (2016). BiGG Models: A platform for integrating, standardizing and sharing genome-scale models. Nucleic acids research, 44(D1), D515-D522.





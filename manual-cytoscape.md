# Diffany Cytoscape plugin ####
## Installation ####
 - If not installed yet, Cytoscape 3 can be downloaded at the [Cytoscape homepage][1]
 - Download the Diffany app file (diffany-0.0.1-RC.jar) at the bottom of the [supplementary data][2]
 - Launch Cytoscape 3
 - In the Cytoscape menu: Apps -> App Manager
 - In the App Manager: Install from File... 
 - Locate the jar file you downloaded above
 - If installed correctly, the "Currently Installed" tab of the App Manager should now display "be.svlandeg.diffany" as "Installed"
 
[1]: http://cytoscape.org/
[2]: http://bioinformatics.psb.ugent.be/supplementary_data/solan/diffany/

## Overview of the functionality ####
The **Diffany tab** in the side panel allows you to select a network collection which you want to use for the inference of differential networks. Once a specific collection is selected from the drop-down list, you need to configure the exact settings of the Diffany run. A network can be included in the project by selecting the checkbox next to its name. Upon inclusion, the visual style "Diffany - Source" is immediately applied on its view. The correct Reference network can further be specified as well, and the user can chose to generate only differential networks, only consensus networks, or both simultaneously.

Using the pairwise **comparison mode**, the reference network will be compared to all other input networks seperately (1 or more), generating condition-specific differential and consensus networks. By contrast, the one-to-all mode further allows the simultaneous comparison of a reference network against many different conditions, generating an overview of the shared stress response across all these conditions. When there are more than two input networks, the minimum number of required supporting networks can further be defined through a **required networks slider**: selecting the number of all input networks (as by default) will result in the most strict comparison, requiring all conditions to agree on a certain response for it to be included in the differential network. By lowering this threshold, it becomes possible to generate differential networks in which for instance 3 out of 4 conditions share a response, thus generating a more general overview while also allowing for some noise in the input data. Finally, a **weight cutoff** value above 0 can be defined to prevent edges with a too small weight to be included in the generated differential and/or consensus networks.

Clicking the **"Start" button** at the bottom of the tab runs the algorithm. If this button is grey, the project is not well configured. After the run, the generated networks will be shown within Cytoscape, and a popup displays the log file of the run. Note that the differential networks will automatically obtain the visual style "Diffany - Differential".

The Diffany algorithms should work on a plethora of networks, taking into account a few caveats:
 - The networks to be compared should be created within the **same network collection**.
 - Weights should be floating point values specified in the **column (node attribute) "weight"**
 - Explicitly negating edges ("does NOT bind") can be obtained by adding a **column "negated"**, containing boolean values
 - **Equivalent nodes** should have the same name accross the different networks in order to compare their edges.
 - Cytoscape considers two edges with the same id to be the exact same edge, which causes attributes (like weight) to be synced accross networks. Take this into account if this is not the intended behaviour.
 - Interaction types are defined by the default **column "interaction"**. When an interaction type not available in the default edge ontology is chosen, it will be transparently added to the ontology.
However, if you want to exploit the existing edge ontology hierarchy, use any of these **built-in AF types**:  genetic\_interaction, positive\_genetic\_interaction, negative\_genetic\_interaction, synthetic\_lethality, regulation, positive\_regulation, catalysis, negative_regulation, inhibition
 or **PR types**: ppi, colocalization, coexpression, protein\_dna\_binding, transcription, ptm, phosphorylation, dephosphorylation, glycosylation, deglycosylation, acetylation, deacetylation, hydroxylation, dehydroxylation, ubiquitination, deubiquitination, methylation, demethylation

## Example networks ####
There are two **examples** included in the app. These can be loaded from the menu: Apps -> Diffany -> Examples. The first one contains two small artificial input networks inspired by the review of Ideker et al., 2012. The second example illustrates the application of Diffany to the study on plant abiotic stress, depicting the reference network (no stress) as well as osmotic stress which was induced 1.5h, 3h, 12h and 24h after transfer to a mannitol-containing medium. Selecting one of the examples will create a new Cytoscape collection with the example networks in it.




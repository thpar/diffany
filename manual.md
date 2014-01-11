# Diffany ####
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

## Example networks ####
There are two small examples included in the app. These can be loaded from the menu: Apps -> Diffany -> Examples. 
Clicking one of the examples will create a new Cytoscape collection with the example networks in it.

The Diffany tab in the side panel will now show the newly added collection(s). Select one of the collections to get a list of its networks.
A network gets included in the project by clicking the checkbox next to it. Upon inclusion, the visual style "Diffany - Source" is immediately applied on its view.
One of the networks needs to be selected as reference network (in the examples resp. "Untreated Network" and "Reference Netwerk"). 
Neither of the examples need further options: they both use a cutoff of 0 and the pairwise comparison mode.

Clicking the "Start" button at the bottom of the tab runs the algorithm.
Two new networks should show up: a differential network (with the "Diffany - Differential" style applied) and its counterpart overlap network.
A popup shows the log file of this run.

## User created networks ####
The Diffany algorithms should work on a plethora of networks, taking into account a few caveats:
 - The networks to be compared should be created within the same network collection.
 - Weights should be floating point values in the column "weight"
 - Explicitly negating edges ("does NOT bind") can be done by adding a column "negated", containing boolean values
 - Equivalent nodes should have the same name accross the different networks in order to compare their edges.
 - Cytoscape sees two edges with the same id as the exact same edge, which causes attributes (like weight) to be synced accross networks. 
Take this into account if this is not the intended behaviour.
 - When there is more than 1 condition-dependent network created, the comparison mode can be set either to pairwise (1-against-1) or 1-against-all.
 - Interaction types are defined by the default column "interaction". When an interaction type not available in the default edge ontology is chosen, it will be transparently added to the ontology.
However, if you want to exploit the existing edge ontology hierarchy, use any of these built-in AF types:  genetic_interaction, positive_genetic_interaction, negative_genetic_interaction, synthetic_lethality, regulation, positive_regulation, catalysis, negative_regulation, inhibition
 or PR types: ppi, colocalization, coexpression, protein_dna_binding, transcription, ptm, phosphorylation, dephosphorylation, glycosylation, deglycosylation, acetylation, deacetylation, hydroxylation, dehydroxylation, ubiquitination, deubiquitination, methylation, demethylation

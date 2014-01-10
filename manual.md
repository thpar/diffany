# Diffany ####
## Installation ####
 - If not installed yet, Cytoscape 3 can be found at http://cytoscape.org/
 - Download the Diffany app file: diffany-0.0.1-RC.jar
 - Launch Cytoscape 3
 - In the Cytoscape menu: Apps -> App Manager
 - In the App Manager: Install from File... 
 - Locate the jar file you downloaded above
 - If installed correctly, the "Currently Installed" tab of the App Manager should now display "be.svlandeg.diffany" as "Installed"
 
## Example networks ####
There are two small examples included in the app. These can be loaded from the menu: Apps -> Diffany -> Examples
Clicking one of the examples will create a new Cytoscape collection with the example networks in it.

The Diffany tab in the side panel will now show the newly added collection(s). A network gets included in the project by clicking the checkbox next to it. Upon inclusion, the visual style "Diffany Source" is immediately applied.
One of the networks needs to be selected as reference network (in the examples resp. "Untreated Network" and "Reference Netwerk"). Neiter of the examples need further options. They both use a cutoff 0 and the Pairwise comparison mode.

Clicking the "Start" button at the bottom of the tab runs the algorithm. Two new networks should show up: a differential network (with the "Diffany Differantial" style applied) and an overlap network.
A popup shows the log file of this run.

## User created networks ####
The Diffany algorithm should work on a plethora of networks, taking into account a few caveats:
 - Weights should be floating point values in the column "weight"
 - Explicitly negating an edge ("does NOT bind") can be done by adding a column "negated", containing boolean values
 - Nodes should have the same name accross the different networks in order to compared their edges.
 - Cytoscape sees two edges with the same id as the exact same edge, which causes attributes (like weight) to be synced accross networks. Take this into account if this is not the intended behaviour.


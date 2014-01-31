# Diffany commandline interface (cli) ####
## Installation ####
 - If not installed yet, java 6 or higher should be installed 
 - Download the Diffany library jar file and add it to your classpath
 - Download the Diffany cli jar file
 - The cli jar is intented for generating text-readable differential networks

## Input data ####
 - The reference network and condition-dependent network both need to be defined in a separate directory with .tab files:
   + **network.tab**: contains a line "Name \t XXX" with XXX the name of the network, and a similar line for the Type: ReferenceNetwork of ConditionNetwork
   + **nodes.tab**: contains one node per line, specifying the node name
   + **edges.tab**: contains one edge per line, in a tab-delimited format: source\_name target\_name  interaction\_type symmetrical/directed  weight  negated/not negated

## Usage parameters ####
 

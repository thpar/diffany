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
   + **edges.tab**: contains one edge per line, in a tab-delimited format: *source\_name* \t *target\_name* \t *interaction\_type* \t *symmetrical/directed* \t *weight* \t *negated/not negated*
   + **conditions.tab**: only needed for the ConditionNetwork. Specifies one condition per line, consisting of its description and a list of tab-delimited ontology terms
   
## Output data ####
 - In a similar format as described above, the differential and the overlapping networks will each be written to a separate output directory with the same .tab files.

## Usage parameters ####
 - Run the jar with **-v** or **--version** (only!) to obtain the version number of the Diffany tool
 - Run the jar with **-h** or **--help** (only!) to obtain this help file:

 - Usage: java -jar Diffany_0.0.1.jar -cond \<dir\> [-conf <arg>] -diff \<dir\> [-l] [-name <arg>] -overlap \<dir\> -ref \<dir\>
   + **-cond**,--conditionsDirectory \<dir\> : the input directory containing the condition-specific network
   + **-conf**,--confidenceMin \<arg\> : the minimum confidence threshold for differential and overlap edges
   + **-diff**,--differentialDirectory \<dir\> : the output directory which will contain the generated differential network
   + **-l**,--log : display a log file after running the algorithm
   + **-name**,--networkOutputName \<arg\> : the name of the generated differential network
   + **-overlap**,--overlappingDirectory \<dir\> : the output directory which will contain the generated overlap network
   + **-ref**,--referenceDirectory \<dir\> : the input directory containing the reference network
 
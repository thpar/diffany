# Diffany commandline interface (cli) ####
## Installation ####
 - If not installed yet, java 6 or higher should be installed 
 - Download the Diffany library jar file and add it to your classpath
 - Download the Diffany cli jar file
 - The cli jar is intented for generating text-readable differential networks

## Input data ####
 - Within one root input directory, the reference network and condition-dependent network(s) need to be defined in separate subdirectories each, containing these .txt files:
   + **network.txt**: contains the ID, name and type of the network, as well as its pre-defined node attributes
   + **nodes.txt**: contains one node per line, specifying the node ID and symbol, as well as the values for the node attributes
   + **edges.txt**: contains one edge per line, in a tab-delimited format: *source\_ID* \t *target\_ID* \t *interaction\_type* \t *symmetrical/directed* \t *weight* \t *negated/affirmative*
   + **conditions.txt**: only needed for the ConditionNetwork. Specifies one condition per line, consisting of its description and an (optional) list of tab-delimited ontology terms
   
## Output data ####
 - In a similar format as described above, the differential and the overlapping networks will each be written to a separate output directory with the same .txt files.

## Usage parameters ####
 - Basic functionality: java -jar Diffany_CL_1.0.0.jar -i \<inputdir\> -o \<outputdir\>

 - Run the jar with **-v** or **--version** (only!) to obtain the version number of the Diffany tool
 - Run the jar with **-h** or **--help** (only!) to obtain the help file with more detailed instructions:
 
 - usage: java -jar Diffany_CL_1.0.0.jar [-c \<arg\>] [-consName \<arg\>] [-consNet \<arg\>] [-diffName \<arg\>] [-diffNet \<arg\>] [-h \<arg\>] -i \<dir\> [-ID \<arg\>] [-l] [-m \<arg\>] -o \<dir\> [-oper \<arg\>]
 + **-c**,--confidence \<arg\>                   the minimum confidence threshold for output edges, as an integer or double (default=0.0)
 + -**consName**,--consensusName \<arg\>         the name of the generated consensus network
 + -**consNet**,--consensusNetworks \<arg\>      whether or not to calculate consensus networks: yes or no (default=yes)
 + -**diffName**,--differentialName \<arg\>      the name of the generated differential network
 + -**diffNet**,--differentialNetworks \<arg\>   whether or not to calculate differential networks: yes or no (default=yes)
 + -**h**,--skipHeader \<arg\>                   whether or not to skip the first line (header) in the network .txt files (default=yes)
 + -**i**,--inputDir \<dir\>                     the input directory containing the reference and condition-specific networks
 + -**ID**,--outputID \<arg\>                    the first ID that will be used for the generated networks
 + -**l**,--log                                display a progress/log file during the run
 + -**m**,--mode \<arg\>                         the mode of comparison: pairwise or all (default=all)
 + -**o**,--outputDir \<dir\>                    the output directory which will contain the generated differential/consensus networks
 + -**oper**,--operator \<arg\>                  the operator used to create consensus edges: min or max (default=min)
 

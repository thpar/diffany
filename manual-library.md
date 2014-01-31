# Diffany library ####
## Installation ####
 - If not installed yet, java 6 or higher should be installed 
 - Download the Diffany library jar file and add it to your classpath
 - The library jar is intented to be used as an external to other custom projects

## Networks ####
 
 - cf. package be.svlandeg.diffany.concepts
 - Within Diffany, there are 4 types of networks, each inheriting from the Network class. They are defined by a name, set of Edges and a set of Nodes. Additionally, a NodeMapper object is required to define when two Nodes within or across networks can be considered to be equal
    -- ReferenceNetwork: a 'static' input network that is used as a reference, i.e. the interactome under unspecified/unknown/wild-type conditions
    -- ConditionNetwork: a 'static' input network that represents an interactome under a specific (set of) Condition(s). A Condition is defined by a description and (optionally) a set of ontology terms.

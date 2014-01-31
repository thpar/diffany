# Diffany library ####
## Installation ####
 - If not installed yet, java 6 or higher should be installed 
 - Download the Diffany library jar file and add it to your classpath
 - The library jar is intented to be used as an external to other custom projects

## Networks ####
 
 - cf. package be.svlandeg.diffany.concepts
 - Within Diffany, there are 4 types of networks, each inheriting from the Network class. They are defined by a name, set of Edges and a set of Nodes. Additionally, a NodeMapper object is required to define when two Nodes within the network can be considered to be equal
    + ReferenceNetwork: a 'static' input network that is used as a reference, i.e. the interactome under unspecified/unknown/wild-type conditions
    + ConditionNetwork: a 'static' input network that represents an interactome under a specific (set of) Condition(s). A Condition is defined by a description and (optionally) a set of ontology terms.
    + DifferentialNetwork: the 'output' network that summarizes the rewiring events from the reference network to the condition-specific network(s)
    + OverlappingNetwork: the 'output' network that summarizes the overlap between a set of original networks, such as the reference network and the condition-specific network(s). Can be accessed as a field of DifferentialNetwork, and can be seen as its logical counterpart.
    

## Semantics ####
 - cf. package be.svlandeg.diffany.semantics
 - A NodeMapper needs to be defined to define when two Nodes across networks are equal. By default, this is determined by comparing the canonical node names through Node.getName(true) in DefaultNodeMapper. This can be customized by making a custom implementation of NodeMapper.
 - An EdgeOntology provides the semantic interpretation of the interaction (edge) types. Through this ontology, it becomes possible to define synonyms/abbreviations (e.g. 'post-translational modification' and 'ptm') as well as sub- and superclasses (e.g. 'phosphorylation' and 'ptm'). The EdgeOntology heavily defines the output of the differential algorithms. A comprehensive set of interaction types is already defined in DefaultEdgeOntolgy and can be further extended.

## Algorithm ####

 - cf. package be.svlandeg.diffany.algorithms
 - The class CalculateDiff provides the methods to generate the differential networks from the input networks. Optional arguments are the minimal weight threshold (default 0.0) and the name of the output networks (by default 'diff\_XXX' and 'overlap\_XXX'
 - As part of the CalculateDiff algorithms, the methods from NetworkCleaning and Unification will be called.

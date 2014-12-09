# Diffany library ####

This readme file explains the basics for installation and usage of the Diffany library code. Not all packages and functionality is explained, but much more information can be found in the javadoc annotations of the original source files.

## Installation ####
 - If not installed yet, java 6 or higher should be installed 
 - Download the Diffany library jar file and add it to your classpath
 - The library jar is intented to be used as an external library to other custom projects

## Networks ####
 - cf. package *be.svlandeg.diffany.core.networks*
 - Within Diffany, there are 4 types of networks, each inheriting from the Network class. They are defined by a name, set of Edges and a set of Nodes. 
    + **ReferenceNetwork**: a 'static' input network that is used as a reference, i.e. the interactome under unspecified/unknown/wild-type conditions.
    + **ConditionNetwork**: a 'static' input network that represents an interactome under a specific (set of) Condition(s). A Condition is defined by a description and (optionally) a set of ontology terms.
    + **DifferentialNetwork**: the 'output' network that summarizes the rewiring events from the reference network to the condition-specific network(s)
    + **OverlappingNetwork**: the 'output' network that summarizes the overlap between a set of original networks, such as the reference network and the condition-specific network(s). Can be intuitively seen as the logical counterpart of a differential network, though their inference of done independently.
    

## Semantics ####
 - cf. package *be.svlandeg.diffany.core.semantics*
 - A **NodeMapper** facilitates the comparison of two Nodes across networks to determine equality. Currently, nodes with the same unique ID are considered equal, but in the future more complex N-M mappings could be implemented.
 - An **EdgeOntology** provides the semantic interpretation of the interaction (edge) types. Through this ontology, it becomes possible to define synonyms/abbreviations (e.g. 'post-translational modification' and 'ptm') as well as sub- and superclasses (e.g. 'phosphorylation' and 'ptm'). The EdgeOntology heavily defines the output of the differential algorithms. A comprehensive set of interaction types is already defined in DefaultEdgeOntolgy and can be further extended.

## Project ####
 - cf. package *be.svlandeg.diffany.core.project*
 - a **Project** keeps track of all ontologies as well as a set of RunConfigurations, i.e. meaningful combinations of reference and condition-specific networks that can be used as input for the Diffany algorithms.
 - Each **RunConfiguration** should be added to the project, and the unique ID assigned to it can be used to run the configuration and retrieve the log file afterwards.
 - A **Logger** object records all relevant messages during the execution of the network algorithms.

## Algorithm ####
 - cf. package *be.svlandeg.diffany.core.algorithms*
 - The class **CalculateDiff** provides the methods to generate the differential networks from the input networks. Optional arguments are for instance the minimal weight threshold (default 0.0) and the name of the output networks (by default 'diff\_XXX' and 'consensus\_XXX'
 - As part of the CalculateDiff algorithms, the methods from **NetworkCleaning** and **Unification** will be called.
 
## Example code ####

Parameters: String refLocation, String condLocation, String diffLocation, String consensusLocation
	
	
		/** DEFINE THE ONTOLOGIES AND THE PROJECT **/
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project("testProject", eo);

		/** READ THE INPUT NETWORKS **/
		boolean skipHeader = true;
		boolean throwInvalidDirException = true;
		File refDir = new File(refLocation);
		ReferenceNetwork refNet = NetworkIO.readReferenceNetworkFromDir(refDir, skipHeader, throwInvalidDirException);

		File condDir = new File(condLocation);
		ConditionNetwork condNet = NetworkIO.readConditionNetworkFromDir(condDir, skipHeader, throwInvalidDirException);

		/** DEFINE THE RUN PARAMETERS **/
		double cutoff = 0.0;
		boolean cleanInput = true;
		int runID = p.addRunConfiguration(refNet, condNet, cleanInput, null);
		
		/** THE ACTUAL ALGORITHM **/
		CalculateDiff diffAlgo = new CalculateDiff();
		diffAlgo.calculateOneDifferentialNetwork(p, runID, cutoff, null, null, 342, 666, true, null);

		// In this case, there will be exactly one DifferentialNetwork
		RunOutput output = p.getOutput(runID);
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork diffNet = pair.getDifferentialNetwork();
		ConsensusNetwork consensusNet = pair.getConsensusNetwork();

		/** WRITE NETWORK OUTPUT **/
		boolean writeHeaders = true;
		
		File diffDir = new File(diffLocation);
		NetworkIO.writeNetworkToDir(diffNet, diffDir, writeHeaders);

		File consensusDir = new File(consensusLocation);
		NetworkIO.writeNetworkToDir(consensusNet, consensusDir, writeHeaders);

		/** WRITE LOG OUTPUT **/
		Logger logger = p.getLogger(runID);
		for (LogEntry msg : logger.getAllLogMessages())
		{
			System.out.println(msg);
		}
	

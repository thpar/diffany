package be.svlandeg.diffany.core.project;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.algorithms.NetworkCleaning;
import be.svlandeg.diffany.core.algorithms.Unification;
import be.svlandeg.diffany.core.io.ProjectIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/**
 * A project consists of a number of Diffany runs and the ontology settings within one user session.
 * 
 * Specifically, it contains one or more {@link RunConfiguration}s which is a subset of networks that can together
 * be used as input for the Diffany algorithms. 
 * Further, for each ID corresponding to such configuration, there is a RunOutput and a Logger object.
 * 
 * Additionally, a project links to an {@link EdgeOntology} that defines the semantics of edge types,
 * and a {@link NodeMapper} that establishes equality of nodes across networks.
 * 
 * Project data can be saved and loaded through the {@link ProjectIO} class.
 * 
 * @author Sofie Van Landeghem
 */
public class Project
{
	
	protected String name;

	protected TreeEdgeOntology edgeOntology;
	protected NodeMapper nodeMapper;
	
	// configurations by Run IDs
	protected Map<Integer, RunConfiguration> configurations;
	
	// outputs by Run IDs
	protected Map<Integer, RunOutput> runOutputs;
	
	// types by Run IDs
	protected Map<Integer, Boolean> runTypes;	// boolean: can do differential
	
	// loggers by Run ID
	protected Map<Integer, Logger> runLogs;
	
	/**
	 * Create a new project with a default node mapper and a default edge ontology that can interpret the differential edges.
	 * 
	 * @param name the name of this project (not null!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public Project(String name) throws IllegalArgumentException
	{
		this(name, new DefaultEdgeOntology(), new DefaultNodeMapper());
	}

	/**
	 * Create a new project with a node mapper and an edge ontology that can interpret the differential edges.
	 * 
	 * @param name the name of this project (not null!)
	 * @param edgeOntology the edge ontology (not null!)
	 * @param nodeMapper the node mapper (not null!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public Project(String name, TreeEdgeOntology edgeOntology, NodeMapper nodeMapper) throws IllegalArgumentException
	{
		if (name == null)
		{
			String errormsg = "The name of a project should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.name = name;
		
		setEdgeOntology(edgeOntology);
		setNodeMapper(nodeMapper);
		
		configurations = new HashMap<Integer, RunConfiguration>();
		runOutputs = new HashMap<Integer, RunOutput>();
		runTypes = new HashMap<Integer, Boolean>();
		runLogs = new HashMap<Integer, Logger>();
	}
	
	/**
	 * Retrieve all run IDs stored in this project
	 * @return all run IDs
	 */
	public Set<Integer> getAllRunIDs()
	{
		return configurations.keySet();
	}
	
	/**
	 * Get a specific RunConfiguration by its unique run ID
	 * 
	 * @param runID the (unique) run ID of the configuration
	 * @return the corresponding RunConfiguration
	 * @throws IllegalArgumentException if the run ID is invalid (i.e. non-existing)
	 */
	public RunConfiguration getRunConfiguration(int runID) throws IllegalArgumentException
	{
		if (! configurations.containsKey(runID))
		{
			String errormsg = "Unknown run ID " + runID + " in Project " + name;
			throw new IllegalArgumentException(errormsg);
		}
		return configurations.get(runID);
	}
	
	
	/**
	 * Add a new RunConfiguration to this project. 
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * 
	 * @param reference the reference network (not null!)
	 * @param condition the condition-specific networks (not null!)
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 */
	public int addRunConfiguration(ReferenceNetwork reference, ConditionNetwork condition)
	{
		Set<ConditionNetwork> cs = new HashSet<ConditionNetwork>();
		cs.add(condition);
		return addRunConfiguration(reference, cs, 1);
	}
	
	/**
	 * Add a new RunConfiguration to this project, automatically registering the networks to this project.
	 * This configuration will allow the calculation of both differential and overlapping networks.
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * This option will use the default requirement of all networks having an overlapping edge before it can be include in the OverlapNetwork.
	 * 
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 */
	public int addRunConfiguration(ReferenceNetwork reference, Set<ConditionNetwork> conditions)
	{
		return addRunConfiguration(reference, conditions, conditions.size() + 1);
	}
	
	/**
	 * Add a new RunConfiguration to this project, automatically registering the networks to this project.
	 * This configuration will allow the calculation of both differential and overlapping networks.
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * 
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * @param overlapNo_cutoff the number of networks that should at least match for overlap to be defined: min. 2, max conditions.size + 1.
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 */
	public int addRunConfiguration(ReferenceNetwork reference, Set<ConditionNetwork> conditions, int overlapNo_cutoff)
	{
		Logger logger = new Logger();
		logger.log("Analysing the reference and condition-specific network(s) ");
		
		registerSourceNetwork(reference, logger);
		ReferenceNetwork cleanRef = new NetworkCleaning(logger).fullInputRefCleaning(reference, nodeMapper, edgeOntology);
		
		Set<ConditionNetwork> cleanConditions = new HashSet<ConditionNetwork>();
		for (ConditionNetwork conNet : conditions)
		{
			registerSourceNetwork(conNet, logger);
			ConditionNetwork cleanCon = new NetworkCleaning(logger).fullInputConditionCleaning(conNet, nodeMapper, edgeOntology);
			cleanConditions.add(cleanCon);
		}
		
		RunConfiguration rc = new RunDiffConfiguration(cleanRef, cleanConditions);
		int nextID;
		if (configurations.keySet().isEmpty())
		{
			nextID = 1;
		}
		else
		{
			nextID = Collections.max(configurations.keySet()) + 1;
		}
		configurations.put(nextID, rc);
		runTypes.put(nextID, true);
		runLogs.put(nextID, logger);
		runOutputs.put(nextID, new RunOutput());
		return nextID;
	}
	
	/**
	 * Add a new RunConfiguration to this project, automatically registering the networks to this project.
	 * This configuration will NOT be able to calculate differential networks; only overlapping networks.
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * This option will use the default requirement of all networks having an overlapping edge before it can be include in the OverlapNetwork.
	 * 
	 * @param inputNetworks all the input networks (at least 1!)
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 */
	public int addRunConfiguration(Set<InputNetwork> networks)
	{
		return addRunConfiguration(networks, networks.size());
	}
	
	/**
	 * Add a new RunConfiguration to this project, automatically registering the networks to this project.
	 * This configuration will NOT be able to calculate differential networks; only overlapping networks.
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * 
	 * @param inputNetworks all the input networks (at least 1!)
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 */
	public int addRunConfiguration(Set<InputNetwork> networks, int overlapNo_cutoff)
	{
		Logger logger = new Logger();
		logger.log("Analysing the input networks ");
		
		Set<InputNetwork> cleanNetworks = new HashSet<InputNetwork>();
		for (InputNetwork inputNet : networks)
		{
			registerSourceNetwork(inputNet, logger);
			InputNetwork cleanNet = new NetworkCleaning(logger).fullInputCleaning(inputNet, nodeMapper, edgeOntology);
			cleanNetworks.add(cleanNet);
		}
		
		RunConfiguration rc = new RunConfiguration(cleanNetworks, overlapNo_cutoff);
		int nextID;
		if (configurations.keySet().isEmpty())
		{
			nextID = 1;
		}
		else
		{
			nextID = Collections.max(configurations.keySet()) + 1;
		}
		configurations.put(nextID, rc);
		runTypes.put(nextID, false);
		runLogs.put(nextID, logger);
		runOutputs.put(nextID, new RunOutput());
		return nextID;
	}
	
	/**
	 * Get the output for a specific run ID. The RunOutput object will be empty if the corresponding run configuration was not deployed.
	 * @return the output of a specific run.
	 */
	public RunOutput getOutput(int runID)
	{
		return runOutputs.get(runID);
	}
	
	/**
	 * Get the logger for a specific run ID. The logs will be empty if the corresponding run configuration was not deployed.
	 * @return the logs of a specific run.
	 */
	public Logger getLogger(int runID)
	{
		return runLogs.get(runID);
	}
	
	/**
	 * Get the type of the run configuration by ID: true if it can calculate differential networks, false otherwise (only overlappingà.
	 * @return the type of a specific run configuration.
	 */
	public boolean isDiffType(int runID)
	{
		return runTypes.get(runID);
	}
	
	
	/**
	 * Register a new source network to this project, updating the EdgeOntology and NodeMapper instances.
	 * A source network is a reference, condition-specific, or an overlapping network.
	 * There is NO check whether or not this network was added previously.
	 * 
	 * @param source the new Network that will be used within this project
	 */
	public void registerSourceNetwork(Network source, Logger logger)
	{
		new Unification(logger).expandEdgeOntology(source.getEdges(), edgeOntology);
	}

	/**
	 * Return the name of this project
	 * @return the name of this project
	 */
	public String getName()
	{
		return name;
	}


	/**
	 * Set the edge ontology for this project.
	 * (currently not a public method - changes to it would influence the differential networks (TODO v3.0))
	 * 
	 * @param edgeOntology the edge ontology (not null!)
	 * @throws IllegalArgumentException if the edgeOntology is null
	 */
	private void setEdgeOntology(TreeEdgeOntology edgeOntology) throws IllegalArgumentException
	{
		if (edgeOntology == null)
		{
			String errormsg = "The edge ontology should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.edgeOntology = edgeOntology;

	}
	
	/**
	 * Get the edge ontology of this project, which can translate edge types to categories and assign semantics to the categories.
	 * 
	 * @return the edge ontology used in this project (should not be null)
	 */
	public TreeEdgeOntology getEdgeOntology()
	{
		return edgeOntology;
	}

	/**
	 * Set the node mapper for this project.
	 * (currently not a public method - changes to it would influence the differential networks (TODO v3.0))
	 * 
	 * @param nodeMapper the node mapper (not null!)
	 * @throws IllegalArgumentException if the nodeMapper is null
	 */
	private void setNodeMapper(NodeMapper nodeMapper) throws IllegalArgumentException
	{
		if (nodeMapper == null)
		{
			String errormsg = "The node mapper should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.nodeMapper = nodeMapper;
	}

	/**
	 * Get the node mapper of this project, which defines equality between nodes of the different networks.
	 * @return the node mapper of this project (should not be null)
	 */
	public NodeMapper getNodeMapper()
	{
		return nodeMapper;
	}


}

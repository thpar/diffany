package be.svlandeg.diffany.concepts;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.algorithms.NetworkCleaning;
import be.svlandeg.diffany.io.ProjectIO;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * A project consists of a number of Diffany runs and the ontology settings within one user session.
 * 
 * Specifically, it contains one or more {@link RunConfiguration}s which is a subset of networks that can together
 * be used as input for the Diffany algorithms.
 * Additionally, a project links to an {@link EdgeOntology} that defines the semantics of edge types,
 * and a {@link NodeMapper} that establishes equality of nodes across networks.
 * Finally, a Logger instance keeps a logfile of all runs in the project.
 * 
 * Project data can be saved and loaded through the {@link ProjectIO} class.
 * 
 * @author Sofie Van Landeghem
 */
public class Project
{
	
	protected String name;

	protected EdgeOntology edgeOntology;
	protected NodeMapper nodeMapper;
	
	// runs by Configuration IDs
	protected Map<Integer, RunConfiguration> runs;
	
	// loggers by Configuration ID
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
	public Project(String name, EdgeOntology edgeOntology, NodeMapper nodeMapper) throws IllegalArgumentException
	{
		if (name == null)
		{
			String errormsg = "The name of a project should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.name = name;
		
		setEdgeOntology(edgeOntology);
		setNodeMapper(nodeMapper);
		
		runs = new HashMap<Integer, RunConfiguration>();
		runLogs = new HashMap<Integer, Logger>();
	}
	
	/**
	 * Get a specific RunConfiguration by its unique ID
	 * 
	 * @param configurationID the (unique) ID of the configuration
	 * @return the corresponding RunConfiguration
	 * @throws IllegalArgumentException if the configuration ID is invalid (i.e. non-existing)
	 */
	public RunConfiguration getRunConfiguration(int configurationID) throws IllegalArgumentException
	{
		if (! runs.containsKey(configurationID))
		{
			String errormsg = "Unknown configuration ID " + configurationID + " in Project " + name;
			throw new IllegalArgumentException(errormsg);
		}
		return runs.get(configurationID);
	}
	
	/**
	 * Add a new RunConfiguration to this project. 
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * 
	 * @param reference the reference network (not null!)
	 * @param condition the condition-specific networks (not null!)
	 * 
	 * @return the unique ID assigned to the RunConfiguration in this project
	 */
	public int addRunConfiguration(ReferenceNetwork reference, ConditionNetwork condition)
	{
		Set<ConditionNetwork> cs = new HashSet<ConditionNetwork>();
		cs.add(condition);
		return addRunConfiguration(reference, cs);
	}
	
	/**
	 * Add a new RunConfiguration to this project. 
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * 
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * 
	 * @return the unique ID assigned to the RunConfiguration in this project
	 */
	public int addRunConfiguration(ReferenceNetwork reference, Set<ConditionNetwork> conditions)
	{
		Logger logger = new Logger();
		ReferenceNetwork cleanRef = new NetworkCleaning(logger).fullInputRefCleaning(reference, nodeMapper, edgeOntology, true);
		Set<ConditionNetwork> cleanConditions = new HashSet<ConditionNetwork>();
		for (ConditionNetwork conNet : conditions)
		{
			ConditionNetwork cleanCon = new NetworkCleaning(logger).fullInputConditionCleaning(conNet, nodeMapper, edgeOntology, true);
			cleanConditions.add(cleanCon);
		}
		
		RunConfiguration rc = new RunConfiguration(cleanRef, cleanConditions);
		int nextID;
		if (runs.keySet().isEmpty())
		{
			nextID = 1;
		}
		else
		{
			nextID = Collections.max(runs.keySet()) + 1;
		}
		runs.put(nextID, rc);
		runLogs.put(nextID, logger);
		return nextID;
	}
	
	/**
	 * Get the logger for a specific runconfiguration. The logs will be empty if the runconfiguration was not deployed.
	 * @return the logs of a specific run configuration.
	 */
	public Logger getLogger(int configurationID)
	{
		return runLogs.get(configurationID);
	}
	
	
	/**
	 * Register a new source network to this project, updating the EdgeOntology and NodeMapper instances.
	 * A source network is a reference, condition-specific, or an overlapping network.
	 * There is NO check whether or not this network was added previously.
	 * 
	 * @param source the new Network that will be used within this project
	 */
	public void registerSourceNetwork(Network source)
	{
		for (Edge e : source.getEdges())
		{
			String edgeType = e.getType();
			if (! edgeOntology.isDefinedSourceType(edgeType))
			{
				boolean symmetrical = e.isSymmetrical();
				
				// add the new edge type as its own separate (singleton) source category
				edgeOntology.addSourceCategory(edgeType, symmetrical);
				edgeOntology.addSourceCategoryMapping(edgeType, edgeType, false);
			}
		}
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
	 * (currently not a public method - changes to it would influence the differential networks (TODO v2.0))
	 * 
	 * @param edgeOntology the edge ontology (not null!)
	 * @throws IllegalArgumentException if the edgeOntology is null
	 */
	private void setEdgeOntology(EdgeOntology edgeOntology) throws IllegalArgumentException
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
	public EdgeOntology getEdgeOntology()
	{
		return edgeOntology;
	}

	/**
	 * Set the node mapper for this project.
	 * (currently not a public method - changes to it would influence the differential networks (TODO v2.0))
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

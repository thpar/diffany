package be.svlandeg.diffany.concepts;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import be.svlandeg.diffany.io.ProjectIO;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * A project consists of a number of input and output networks within a certain session. 
 * 
 * It contains one or more {@link RunConfiguration}s which is a subset of networks that can together
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
	
	protected Logger logger;

	protected Map<Integer, RunConfiguration> runs;

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
		
		logger = new Logger();
		runs = new HashMap<Integer, RunConfiguration>();
	}
	
	/**
	 * Get the RunConfiguration by its unique ID
	 * 
	 * @param configurationID the (unique) ID of the configuration
	 * @return the corresponding RunConfiguration
	 * @ throws IllegalArgumentException if the configuration ID is invalid (i.e. non-existing)
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
	 * Add a new RunConfiguration to this project
	 * @param rc the new  new RunConfiguration
	 * @return the unique ID assigned to the RunConfiguration in this project
	 */
	public int addRunConfiguration(RunConfiguration rc)
	{
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
		return nextID;
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
	 * Get the edge ontology of this project, which can translate edge types to categories and assign semantics to the categories.
	 * 
	 * @return the edge ontology used in this project (should not be null)
	 */
	public EdgeOntology getEdgeOntology()
	{
		return edgeOntology;
	}

	/**
	 * Get the node mapper of this project, which defines equality between nodes of the different networks.
	 * @return the node mapper of this project (should not be null)
	 */
	public NodeMapper getNodeMapper()
	{
		return nodeMapper;
	}

	
	/**
	 * Get the logger object, that contains all log messages from the last run.
	 * @return the logger object of this project.
	 */
	public Logger getLogger()
	{
		return logger;
	}

}

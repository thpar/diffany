package be.svlandeg.diffany.concepts;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * A project consists of a selection of networks and all analyses performed on these networks within the current session. 
 * It should contain exactly 1 reference network, at least 1 condition-specific network, 
 * and it may contain 1 or more differential networks.
 * 
 * Additionally, a project links to an ontology that defines the semantics of edge types.
 * 
 * Classes that extend the Project class should provide functionality for saving and loading the data.
 * 
 * @author Sofie Van Landeghem
 */
public class Project
{
	
	protected String name;

	protected EdgeOntology edgeOntology;
	protected NodeMapper nodeMapper;
	
	protected Logger logger;

	protected ReferenceNetwork reference;
	protected Set<ConditionNetwork> conditions;
	protected Set<DifferentialNetwork> differentials;

	/**
	 * Create a new project with a reference network, set of condition-specific networks, 
	 * differential networks, a node mapper and an edge ontology that can interpret the differential edges.
	 * 
	 * @param name the name of this project (not null!)
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * @param differentials the differential networks (not null!)
	 * @param edgeOntology the edge ontology (not null!)
	 * @param nodeMapper the node mapper (not null!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public Project(String name, ReferenceNetwork reference, Set<ConditionNetwork> conditions, Set<DifferentialNetwork> differentials, EdgeOntology edgeOntology, NodeMapper nodeMapper) throws IllegalArgumentException
	{
		if (name == null)
		{
			String errormsg = "The name of a project should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.name = name;
		
		setReference(reference);
		setConditions(conditions);
		setDifferentials(differentials);
		setEdgeOntology(edgeOntology);
		setNodeMapper(nodeMapper);
		
		logger = new Logger();
	}

	/**
	 * Create a new project with a reference network, a set of condition-specific networks, 
	 * a node mapper and  and an edge ontology that can interpret the differential edges. 
	 * The set of differential networks will be empty but can be calculated with the CalculateDiff class.
	 * 
	 * @param name the name of this project (not null!)
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * @param edgeOntology the edge ontology (not null!)
	 * @param nodeMapper the node mapper (not null!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public Project(String name, ReferenceNetwork reference, Set<ConditionNetwork> conditions, EdgeOntology edgeOntology, NodeMapper nodeMapper)
	{
		this(name, reference, conditions, new HashSet<DifferentialNetwork>(), edgeOntology, nodeMapper);
	}

	/**
	 * Create a previous project by loading it from a file location
	 * 
	 * @param location the location where this project was saved previously
	 */
	public Project(String location)
	{
		loadFromFile(location);
		logger = new Logger();
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
	 * Save the project data to a specific file location
	 * 
	 * @param fileLocation the location where the project should be saved
	 */
	public void saveProject(String fileLocation)
	{
		//TODO v1.1: implement!
		throw new UnsupportedOperationException("Saving of project not yet implemented");
	}

	/**
	 * Load the project data from a specific file location. 
	 * Make sure all restrictions on number of networks are respected during the load!
	 *
	 * @param fileLocation the location from where the project should be loaded
	 */
	public void loadFromFile(String fileLocation)
	{
		//TODO v1.1: implement!
		throw new UnsupportedOperationException("Loading of project not yet implemented");
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
	 * Set the condition-specific networks in this project.
	 * (currently not a public method - changes to it would influence the differential networks (TODO v2.0))
	 * 
	 * @param conditions the condition-specific networks (not null or empty!)
	 * @throws IllegalArgumentException if the set is null or empty
	 */
	private void setConditions(Set<ConditionNetwork> conditions) throws IllegalArgumentException
	{
		if (conditions == null || conditions.isEmpty())
		{
			String errormsg = "The set of condition-specific networks should not be null or empty!";
			throw new IllegalArgumentException(errormsg);
		}
		this.conditions = conditions;

	}

	/**
	 * Set the reference network of this project.
	 * (currently not a public method - changes to it would influence the differential networks (TODO v2.0))
	 * 
	 * @param reference the reference network
	 * @throws IllegalArgumentException if the network is null
	 */
	private void setReference(ReferenceNetwork reference) throws IllegalArgumentException
	{
		if (reference == null)
		{
			String errormsg = "The reference network should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.reference = reference;

	}

	/**
	 * Initialize the set of differential networks in this project
	 * 
	 * @param differentials the differential networks (not null!)
	 * @throws IllegalArgumentException if the set is null
	 */
	private void setDifferentials(Set<DifferentialNetwork> differentials) throws IllegalArgumentException
	{
		if (differentials == null)
		{
			String errormsg = "The set of differential networks should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.differentials = differentials;

	}

	/**
	 * Add a differential network to this project.
	 * 
	 * @param differential a new differential network
	 */
	public void addDifferential(DifferentialNetwork differential)
	{
		differentials.add(differential);
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
	 * Get the reference network of this project, against which the condition
	 * dependent network(s) will be compared to.
	 * 
	 * @return the reference network in this project (should not be null)
	 */
	public ReferenceNetwork getReferenceNetwork()
	{
		return reference;
	}

	/**
	 * Get the condition-dependent network(s): 1 or many. 
	 * Should there be 0 condition-dependent networks, the complete project is invalid/incomplete.
	 * 
	 * @return the condition-dependent networks in this project (1 or more, never null or empty))
	 */
	public Collection<ConditionNetwork> getConditionNetworks()
	{
		return conditions;
	}

	/**
	 * Get the differential networks in the project: 0, 1 or more
	 * 
	 * @return the differential networks in this project (if any, otherwise empty set, but never null)
	 */
	public Collection<DifferentialNetwork> getDifferentialNetworks()
	{
		return differentials;
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

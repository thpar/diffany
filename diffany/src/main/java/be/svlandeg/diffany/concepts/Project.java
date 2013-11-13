package be.svlandeg.diffany.concepts;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.semantics.EdgeOntology;

/**
 * A project consists of a selection of networks and all analyses performed on
 * these networks within the current session. It should contain exactly 1
 * reference network, at least 1 condition-specific network, and it may contain
 * 1 or more differential networks.
 * 
 * Additionally, a project links to an ontology that defines the semantics of
 * edge types.
 * 
 * Classes that extends the Project class should provide functionality for saving and loading the data.
 * 
 * @author Sofie Van Landeghem
 * 
 */
public abstract class Project
{

	protected EdgeOntology edgeOntology;

	protected ReferenceNetwork reference;
	protected Set<ConditionNetwork> conditions;
	protected Set<DifferentialNetwork> differentials;

	/**
	 * Create a new project with a reference network, set of condition-specific networks, differential networks 
	 * and an edge ontology that can interpret the differential edges. 
	 * 
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * @param differentials the differential networks (not null!)
	 * @param edgeOntology the edge ontology (not null!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public Project(ReferenceNetwork reference, Set<ConditionNetwork> conditions, Set<DifferentialNetwork> differentials, EdgeOntology edgeOntology) throws IllegalArgumentException
	{
		setReference(reference);
		setConditions(conditions);
		setDifferentials(differentials);
		setEdgeOntology(edgeOntology);
	}
	
	/**
	 * Create a new project with a reference network, a set of condition-specific networks, 
	 * and an edge ontology that can interpret the differential edges.
	 * The set of differential networks will be empty.
	 *  
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * @param edgeOntology the edge ontology (not null!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public Project(ReferenceNetwork reference, Set<ConditionNetwork> conditions, EdgeOntology edgeOntology)
	{
		this(reference, conditions, new HashSet<DifferentialNetwork>(), edgeOntology);
	}
	
	/**
	 * Create a previous project by loading it from a file location
	 * @param location the location where this project was saved previously
	 */
	public Project(String location)
	{
		loadFromFile(location);
	}
	
	/**
	 * Save the project data to a specific file location
	 * @param fileLocation the location where the project should be saved
	 */
	public abstract void saveProject(String fileLocation);
	
	/**
	 * Load the project data from a specific file location. Make sure all restrictions on number of networks are respected during the load!
	 * @param fileLocation
	 */
	public abstract void loadFromFile(String fileLocation);

	/**
	 * Set the edge ontology for this project 
	 * (currently not a public method - changes to it would influence the differential networks (TODO)).
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
	 * Set the condition-specific networks in this project
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
	 * Set the reference network
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
	 * Add a differential network to this project
	 * @param differential a new differential network
	 */
	public void addDifferential(DifferentialNetwork differential)
	{
		differentials.add(differential);
	}


	/**
	 * Get the edge ontology of this project, which can translate edge types to
	 * categories and assign semantics to the categories
	 * 
	 * @return the edge ontology used in this project (should not be null)
	 */
	public EdgeOntology getEdgeOntology()
	{
		return edgeOntology;
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
	 * Get the condition-dependent network(s): 1 or many. Should there be 0
	 * condition-dependent networks, the complete project is invalid/incomplete.
	 * 
	 * @return the condition-dependent networks in this project (1 or more,
	 *         never null or empty))
	 */
	public Collection<ConditionNetwork> getConditionNetworks()
	{
		return conditions;
	}

	/**
	 * Get the differential networks in the project: 0, 1 or more
	 * 
	 * @return the differential networks in this project (if any, otherwise
	 *         empty set, but never null)
	 */
	public Collection<DifferentialNetwork> getDifferentialNetworks()
	{
		return differentials;
	}

}

package be.svlandeg.cytoscape.diffany.concepts;

import java.util.HashSet;
import java.util.Set;

/**
 * This class describes an experimental condition that can be linked to a certain network.
 * A condition is described in free text, but may also be associated to ontology terms.
 * 
 * @author Sofie Van Landeghem
 *
 */
public class Condition
{
	
	protected String description;
	protected Set<String> ontologies;
	
	/**
	 * Create a new condition that can be linked to a network.
	 * 
	 * @param description free-text description of the condition
	 * @param ontologies set of corresponding ontology terms
	 */
	public Condition(String description, Set<String> ontologies)
	{
		this.description = description;
		setOntologyTerms(ontologies);
	}
	
	/**
	 * Create a new condition that can be linked to a network.
	 * The set of ontology terms will be initialized to an empty set.
	 * 
	 * @param description free-text description of the condition
	 */
	public Condition(String description)
	{
		this(description, new HashSet<String>());
	}
	
	
	/**
	 * Return a free-text description of the condition
	 * @return the description of the condition
	 */
	public String getDescription()
	{
		return description;
	}
	
	/**
	 * Return the corresponding ontology terms for this condition
	 * @return the ontology terms for this condition (may be empty)
	 */
	public Set<String> getOntologyTerms()
	{
		return ontologies;
	}
	
	/**
	 * Define the set of ontologies corresponding to this condition 
	 * (overwrites all previously defined terms).
	 * @param ontologies the set of ontology terms
	 */
	public void setOntologyTerms(Set<String> ontologies)
	{
		this.ontologies = ontologies;
	}
	
	/**
	 * Add a number of ontology terms to the set of previously defined terms.
	 * @param ontologies the set of ontology terms
	 */
	public void addOntologyTerm(String ontology)
	{
		ontologies.add(ontology);
	}

}

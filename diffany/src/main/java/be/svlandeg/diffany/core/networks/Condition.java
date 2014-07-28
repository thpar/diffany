package be.svlandeg.diffany.core.networks;

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
	 * @param ontologies set of corresponding ontology terms (should not be null)
	 * 
	 * @throws IllegalArgumentException when the description or the set of ontologies is null
	 */
	public Condition(String description, Set<String> ontologies) throws IllegalArgumentException
	{
		if (description == null)
		{
			String errormsg = "The description should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.description = description;
		setOntologyTerms(ontologies);
	}
	
	/**
	 * Create a new condition that can be linked to a network.
	 * The set of ontology terms will be initialized to an empty set.
	 * 
	 * @param description free-text description of the condition
	 * @throws IllegalArgumentException when the description is null
	 */
	public Condition(String description) throws IllegalArgumentException
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
	 * @return the ontology terms for this condition (may be empty but should not be null)
	 */
	public Set<String> getOntologyTerms()
	{
		return ontologies;
	}
	
	/**
	 * Define the set of ontologies corresponding to this condition 
	 * (overwrites all previously defined terms).
	 * 
	 * @param ontologies the set of ontology terms (should not be null)
	 * @throws IllegalArgumentException when the description or the set of ontologies is null
	 */
	public void setOntologyTerms(Set<String> ontologies) throws IllegalArgumentException
	{
		if (ontologies == null)
		{
			String errormsg = "The set of ontologies should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.ontologies = ontologies;
	}
	
	/**
	 * Add a number of ontology terms to the set of previously defined terms.
	 * @param ontology a new ontology term
	 */
	public void addOntologyTerm(String ontology)
	{
		ontologies.add(ontology);
	}

}

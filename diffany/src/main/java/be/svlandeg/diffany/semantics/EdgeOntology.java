package be.svlandeg.diffany.semantics;

import java.util.HashMap;
import java.util.Map;

/**
 * This class takes care of the semantic interpretation of different edge types
 * and their corresponding categories. It can define the category of a
 * differential edge from the categories of the two original edges.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class EdgeOntology
{

	private Map<String, String> mapEdgeToCategory;

	/**
	 * Create an empty ontology, with no edge categories or any mapping defined
	 */
	public EdgeOntology()
	{
		mapEdgeToCategory = new HashMap<String, String>();
	}

	/**
	 * Get the semantic category of a certain edge type. Matching is done independent of upper/lower casing.
	 * 
	 * @param edgeType
	 *            the original type of the edge in a network
	 * @return the semantic category of that edge or null if it is not mapped in this ontology
	 */
	public String getCategory(String edgeType)
	{
		return mapEdgeToCategory.get(edgeType.toLowerCase());
	}

	/**
	 * Create a new mapping from edge type to category. Matching is done independent of upper/lower casing.
	 * 
	 * @param edgeType the original edge type - should not have be defined in this ontology before
	 * @param category the category to be assigned to this edge type
	 * @param overwrite determines whether or not this function may overwrite previous mappings of the same edge type
	 * @throws IllegalArgumentException when the edge type was already mapped in this ontology and overwrite is off
	 */
	public void addCategoryMapping(String edgeType, String category, boolean overwrite) throws IllegalArgumentException
	{
		if (!overwrite && mapEdgeToCategory.containsKey(edgeType.toLowerCase()))
		{
			String errormsg = "The provided edge type is already mapped to a category!";
			throw new IllegalArgumentException(errormsg);
		}
		mapEdgeToCategory.put(edgeType.toLowerCase(), category.toLowerCase());
	}
	

	/**
	 * Method that defines the differential edge category from the corresponding edge categories in the reference and condition-specific networks.
	 * Category names are treated independent of upper/lower casing.
	 * Should only return null when the edge should be deleted (i.e. not present in differential network!)
	 * 
	 * @param referenceCategory the category of the edge in the reference network
	 * @param conditionCategory the category of the edge in the condition-specific network
	 * @return the category of the edge in the differential network, or null when there should be no such edge
	 * @throws IllegalArgumentException when the category of the reference or condition-specific edge does not exist in this ontology
	 */
	public abstract String getDifferentialCategory(String referenceCategory, String conditionCategory) throws IllegalArgumentException;
	

}

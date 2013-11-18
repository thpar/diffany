package be.svlandeg.diffany.semantics;

import java.util.HashMap;
import java.util.Map;

import be.svlandeg.diffany.concepts.Edge;

/**
 * This class takes care of the semantic interpretation of different edge types
 * and their corresponding categories. It can define the category of a
 * differential edge from the categories of the two original edges.
 * 
 * TODO: mapping between categories (is-a).
 * 
 * @author Sofie Van Landeghem
 */
public abstract class EdgeOntology
{

	public static final String VOID_EDGE = "VOID";
	protected Map<String, String> mapEdgeToCategory;
	protected Map<String, Boolean> allCategories;

	/**
	 * Create an empty ontology, with no edge categories or any mapping defined
	 */
	public EdgeOntology()
	{
		removeAllCategoriesAndMappings();
	}

	/**
	 * Get the semantic category of a certain edge type. Matching is done independent of upper/lower casing.
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return the semantic category of that edge or null if it is not mapped in this ontology
	 */
	public String getCategory(String edgeType)
	{
		return mapEdgeToCategory.get(edgeType.toLowerCase());
	}
	
	/**
	 * Determine whether an edge category should be defined as symmetrical or not.
	 * 
	 * @param category the semantic category of an edge
	 * @return whether or not it is symmetrical
	 * @throws IllegalArgumentException if the category is not defined in this ontology
	 */
	public boolean isSymmetrical(String category) throws IllegalArgumentException
	{
		if (! allCategories.containsKey(category))
		{
			String errormsg = "The provided edge category ('" + category + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		return allCategories.get(category);
	}
	

	/**
	 * Create a new mapping from edge type to category. Matching is done independent of upper/lower casing.
	 * 
	 * @param edgeType the original edge type - should not have been defined in this ontology before
	 * @param category the category to be assigned to this edge type
	 * @param overwrite determines whether or not this function may overwrite previous mappings of the same edge type
	 * @throws IllegalArgumentException when the edge type was already mapped in this ontology and overwrite is off,
	 * 	or when the specified category is not part of this ontology
	 */
	public void addCategoryMapping(String edgeType, String category, boolean overwrite) throws IllegalArgumentException
	{
		if (! allCategories.containsKey(category))
		{
			String errormsg = "The provided edge category ('" + category + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!overwrite && mapEdgeToCategory.containsKey(edgeType.toLowerCase()))
		{
			String errormsg = "The provided edge type is already mapped to a category!";
			throw new IllegalArgumentException(errormsg);
		}
		mapEdgeToCategory.put(edgeType.toLowerCase(), category.toLowerCase());
	}
	
	/**
	 * Remove all categories and type-categoryMappings
	 */
	public void removeAllCategoriesAndMappings()
	{
		allCategories = new HashMap<String, Boolean>();
		allCategories.put(VOID_EDGE, true);
		mapEdgeToCategory = new HashMap<String, String>();
	}
	
	/**
	 * Remove all type-category mappings (keeping all defined categories).
	 */
	public void removeAllCategoryMappings()
	{
		mapEdgeToCategory = new HashMap<String, String>();
	}
	

	/**
	 * Method that defines the differential edge category from the corresponding edge categories in the reference and condition-specific networks.
	 * Category names are treated independent of upper/lower casing.
	 * Returns VOID_EDGE when the edge should be deleted (i.e. not present in differential network).
	 * 
	 * @param referenceCategory the category of the edge in the reference network (or VOID_EDGE when non-existing)
	 * @param conditionCategory the category of the edge in the condition-specific network (or VOID_EDGE when non-existing)
	 * @return the category of the edge in the differential network, or VOID_EDGE when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the category of the reference or condition-specific edge does not exist in this ontology
	 */
	public abstract String getDifferentialCategory(String referenceCategory, String conditionCategory) throws IllegalArgumentException;
	
	
}

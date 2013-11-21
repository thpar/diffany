package be.svlandeg.diffany.semantics;

import java.util.*;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This class takes care of the semantic interpretation of different edge types
 * and their corresponding categories. It can define the category of a
 * differential edge from the categories of the two original edges.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class EdgeOntology
{

	public static final String VOID_TYPE = EdgeDefinition.getVoidEdge().getType().toLowerCase();

	private Map<String, String> mapEdgeToCategory;
	private Set<String> allCategories;

	/**
	 * Create an empty ontology, with no edge categories or any mapping defined
	 */
	public EdgeOntology()
	{
		removeAllCategoriesAndMappings();
	}

	/**
	 * Method that defines the differential edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns EdgeDefinition.getVoidEdge() when the edge should be deleted (i.e. not present in differential network).
	 * 
	 * @param referenceEdge the edge definition in the reference network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param conditionEdge the edge definition in the condition-specific network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential network
	 * 
	 * @return the edge definition in the differential network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public abstract EdgeDefinition getDifferentialEdge(EdgeDefinition referenceEdge, EdgeDefinition conditionEdge, double cutoff) throws IllegalArgumentException;

	/**
	 * Method that defines the overlapping edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns EdgeDefinition.getVoidEdge() when the edge should be deleted (i.e. not present in the shared network).
	 * 
	 * @param referenceEdge the edge definition in the reference network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param conditionEdge the edge definition in the condition-specific network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param cutoff the minimal value of a resulting edge for it to be included in the shared network
	 * 
	 * @return the edge definition in the shaerd network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public abstract EdgeDefinition getSharedEdge(EdgeDefinition referenceEdge, EdgeDefinition conditionEdge, double cutoff) throws IllegalArgumentException;

	/**
	 * Get the semantic category of a certain edge type. Matching is done independent of upper/lower casing.
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return the semantic category of that edge or null if it is not mapped in this ontology
	 */
	protected String getCategory(String edgeType)
	{
		return mapEdgeToCategory.get(edgeType.toLowerCase());
	}

	/**
	 * Add a category (casing independent)
	 * @param category the category that should be added to this ontology
	 * @throws IllegalArgumentException when the category is null
	 */
	protected void addCategory(String category) throws IllegalArgumentException
	{
		if (category == null)
		{
			String errormsg = "The category should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		allCategories.add(category.toLowerCase());
	}

	/**
	 * Add a number of categories (casing independent)
	 * @param categories the categories that should be added to this ontology
	 */
	protected void addCategories(Set<String> categories)
	{
		for (String c : categories)
		{
			allCategories.add(c.toLowerCase());
		}
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
		if (!allCategories.contains(category.toLowerCase()))
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
	protected void removeAllCategoriesAndMappings()
	{
		allCategories = new HashSet<String>();
		addCategory(VOID_TYPE);

		mapEdgeToCategory = new HashMap<String, String>();
		addCategoryMapping(VOID_TYPE, VOID_TYPE, false);
	}

	/**
	 * Remove all type-category mappings (keeping all defined categories).
	 */
	protected void removeAllCategoryMappings()
	{
		mapEdgeToCategory = new HashMap<String, String>();
	}

}

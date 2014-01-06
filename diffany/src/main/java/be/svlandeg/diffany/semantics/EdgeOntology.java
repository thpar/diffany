package be.svlandeg.diffany.semantics;

import java.awt.Color;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;
import be.svlandeg.diffany.concepts.VisualEdgeStyle;
import be.svlandeg.diffany.concepts.VisualEdgeStyle.ArrowHead;

/**
 * This class takes care of the semantic interpretation of different edge types and their corresponding categories
 * in the 'source' networks, i.e. the reference and condition-specific networks used as input.
 * It can define a differential edge by inspecting the two original edges (implemented by overriding classes).
 * Additionally, an overlapping edge can be generated by taking into account the edge type-to-category mapping.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class EdgeOntology
{

	public static final String VOID_TYPE = EdgeDefinition.getVoidEdge().getType().toLowerCase();

	protected Map<String, String> mapSourceTypeToCategory;
	
	protected Map<String, Boolean> allSourceCategories;
	protected Set<String> allDiffCategories;
	
	protected Map<String, String> differential_translations;

	/**
	 * Create an empty ontology, with no edge categories or any mapping defined
	 */
	public EdgeOntology()
	{
		removeAllCategoriesAndMappings();
		differential_translations = new HashMap<String, String>();
	}
	
	
	/**
	 * Method that defines the differential edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns EdgeDefinition.getVoidEdge() when the edge should be deleted (i.e. not present in differential network).
	 * 
	 * @param refEdge the edge definition in the reference network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param conEdges the edge definitions in the condition-specific networks (can be EdgeDefinition.getVoidEdge() when non-existing)
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential network
	 * 
	 * @return the edge definition in the differential network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public abstract EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, Set<EdgeDefinition> conEdges, double cutoff) throws IllegalArgumentException;

	/**
	 * Method that defines the overlapping edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns EdgeDefinition.getVoidEdge() when the edge should be deleted (i.e. not present in the overlapping network).
	 * 
	 * @param edges the original edge definitions (can contain EdgeDefinition.getVoidEdge()), should not be empty!
	 * @param cutoff the minimal value of a resulting edge for it to be included in the overlapping network
	 * @param minOperator whether or not to take the minimum of the edge weights - if false, the maximum is taken
	 * 
	 * @return the edge definition in the overlapping network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public abstract EdgeDefinition getOverlapEdge(Set<EdgeDefinition> edges, double cutoff, boolean minOperator) throws IllegalArgumentException;
	
	/**
	 * Define a translation of a type translation into a simplification, 
	 * e.g. "negativeregulation_to_positiveregulation" turns into "increase_regulation"
	 * 
	 * @param transformation the original type change
	 * @param simplification a simplified version of the string, representing the same change
	 * @throws IllegalArgumentException when the transformation was already previously entered
	 */
	public void putDifferentialTranslation(String transformation, String simplification) throws IllegalArgumentException
	{
		if (differential_translations.containsKey(transformation))
		{
			String errormsg = "The provided translation rule ('" + transformation + "') generates a conflict!";
			throw new IllegalArgumentException(errormsg);
		}
		differential_translations.put(transformation, simplification);
	}
	
	/**
	 * Return the translation of a type into its simplification.
	 * @param transformation the original type change
	 * @return its translation (or the same string if there is no mapped translation)
	 */
	public String getDifferentialTranslation(String transformation)
	{
		String simpl = transformation;
		if (differential_translations.containsKey(transformation))
		{
			simpl = differential_translations.get(transformation);
		}
		return simpl;
	}

	/**
	 * Private method that increments a specific type in a category-to-count hashmap.
	 * 
	 * @param catByCount the hashmap holding the counts
	 * @param cat the category that needs to be incremented
	 */
	protected void addOne(Map<String, Integer> catByCount, String cat)
	{
		if (!catByCount.containsKey(cat))
			catByCount.put(cat, 0);
		catByCount.put(cat, catByCount.get(cat) + 1);
	}

	/**
	 * Private method that holds the max depth value of a category mapping:
	 * each time a new depth is found, it influences the new maximum of the map
	 * 
	 * @param catByDepth the hashmap holding the depths
	 * @param depth the newly found depth of the category
	 * @param cat the category that needs to be adjusted in depth
	 */
	protected void recordMaxDepth(Map<String, Integer> catByDepth, String cat, int depth)
	{
		int previousDepth = 0;
		if (catByDepth.containsKey(cat))
			previousDepth = catByDepth.get(cat);
		catByDepth.put(cat, Math.max(depth, previousDepth));
	}

	/**
	 * Get the semantic category of a certain edge type.
	 * Matching is done independent of upper/lower casing.
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return the semantic category of that edge or null if it is not mapped in this ontology
	 */
	public String getSourceCategory(String edgeType)
	{
		return mapSourceTypeToCategory.get(edgeType.toLowerCase());
	}

	/**
	 * Check whether a certain edge type is present in this ontology.
	 * Matching is done independent of upper/lower casing.
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return whether or not this edge type is present in this ontology
	 */
	protected boolean isDefinedSourceType(String edgeType)
	{
		if (edgeType == null)
		{
			return false;
		}
		String cat = getSourceCategory(edgeType);
		return isDefinedSourceCat(cat);
	}
	
	/**
	 * Check whether a certain edge category is present in this ontology.
	 * Matching is done independent of upper/lower casing.
	 * 
	 * @param category the original category of the edge in a network
	 * @return whether or not this edge category is present in this ontology
	 */
	protected boolean isDefinedSourceCat(String category)
	{
		if (category == null)
			 return false;
		return allSourceCategories.containsKey(category);
	}
	
	/**
	 * Check whether a certain edge category is present in this ontology.
	 * Matching is done independent of upper/lower casing.
	 * 
	 * @param category the category of the edge in a (differential) network
	 * @return whether or not this edge category is present in this ontology
	 */
	protected boolean isDefinedDiffCategory(String category)
	{
		if (category == null)
		{
			return false;
		}
		return allDiffCategories.contains(category);
	}
	
	/**
	 * Define the Color object of an edge in a differential network, by edge category.
	 * 
	 * @param category the category of the edge interaction
	 * @return the color of the edge
	 */
	protected abstract Color getDifferentialEdgeColor(String category);
	
	/**
	 * Define the ArrowHead object of an edge in a differential network, by edge category.
	 * 
	 * @param category the category of the edge interaction
	 * @return the arrowhead of the edge
	 */
	protected abstract ArrowHead getDifferentialEdgeArrowHead(String category);
	
	
	/**
	 * Define the full visual style of an edge in a differential network, by edge category.
	 * 
	 * @param category the category of the edge interaction
	 * @return a VisualEdgeStyle object which specifies how the edge should be drawn
	 */
	public VisualEdgeStyle getDifferentialEdgeStyle(String category)
	{
		return new VisualEdgeStyle(getDifferentialEdgeColor(category), getDifferentialEdgeArrowHead(category));
	}
	
	/**
	 * Define the Color object of an edge in a 'normal' network (reference, condition-dependent or overlap).
	 * 
	 * @param edgeType the type of the edge interaction
	 * @return the color of the edge
	 */
	protected abstract Color getSourceEdgeColor(String edgeType);
	
	/**
	 * Define the ArrowHead object of an edge in a 'normal' network (reference, condition-dependent or overlap).
	 * 
	 * @param edgeType the type of the edge interaction
	 * @return the arrowhead of the edge
	 */
	protected abstract ArrowHead getSourceEdgeArrowHead(String edgeType);
	
	/**
	 * Define the visual style of an edge in a 'normal' network (reference, condition-dependent or overlap).
	 * 
	 * @param edgeType the type of the edge interaction
	 * @return a VisualEdgeStyle object which specifies how the edge should be drawn
	 */
	public VisualEdgeStyle getSourceEdgeStyle(String edgeType)
	{
		return new VisualEdgeStyle(getSourceEdgeColor(edgeType), getSourceEdgeArrowHead(edgeType));
	}


	/**
	 * Return all source (input) categories present in this ontology.
	 * @return the set of categories mapped in this ontology.
	 */
	public Set<String> getAllSourceCategories()
	{
		return allSourceCategories.keySet();
	}

	/**
	 * Return all differential categories present in this ontology.
	 * @return the set of categories mapped in this ontology.
	 */
	public Set<String> getAllDiffCategories()
	{
		return allDiffCategories;
	}
	

	/**
	 * Add a differential category (casing independent)
	 * @param category the category that should be added to this ontology
	 * @throws IllegalArgumentException when the category is null
	 */
	protected void addDiffCategory(String category) throws IllegalArgumentException
	{
		if (category == null)
		{
			String errormsg = "The differential category should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		allDiffCategories.add(category.toLowerCase());
	}

	/**
	 * Add a number of differential categories (casing independent)
	 * @param categories the categories that should be added to this ontology
	 */
	protected void addDiffCategories(Set<String> categories)
	{
		for (String c : categories)
		{
			addDiffCategory(c);
		}
	}
	
	/**
	 * Retrieve the symmetry state of a source edgeType in this ontology
	 * @param edgeType the source edgeType
	 * @return the symmetry state of the source category
	 */
	public boolean isSymmetricalSourceType(String edgeType)
	{
		if (edgeType == null)
		{
			String errormsg = "The source category should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (! isDefinedSourceType(edgeType))
		{
			String errormsg = "The source edgeType is not defined in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		String cat = getSourceCategory(edgeType);
		return allSourceCategories.get(cat.toLowerCase());
	}

	/**
	 * Add a source category (casing independent)
	 * @param category the category that should be added to this ontology
	 * @param symmetric whether or not the category is symmetric
	 * @throws IllegalArgumentException when the category is null
	 */
	protected void addSourceCategory(String category, boolean symmetric) throws IllegalArgumentException
	{
		if (category == null)
		{
			String errormsg = "The source category should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (allSourceCategories.containsKey(category))
		{
			String errormsg = "The source category was already previously defined!";
			throw new IllegalArgumentException(errormsg);
		}
		allSourceCategories.put(category.toLowerCase(), symmetric);
	}

	/**
	 * Add a number of source categories (casing independent)
	 * @param categories the categories that should be added to this ontology
	 */
	protected void addSourceCategories(Map<String, Boolean> categories)
	{
		for (String c : categories.keySet())
		{
			addSourceCategory(c, categories.get(c));
		}
	}
	

	/**
	 * Create a new mapping from edge type to category. Matching is done independent of upper/lower casing.
	 * 
	 * @param edgeType the original edge type - should not have been defined in this ontology before
	 * @param category the category to be assigned to this edge type
	 * @param overwrite determines whether or not this function may overwrite previous mappings of the same edge type
	 * @throws IllegalArgumentException when the edge type was already mapped in this ontology and overwrite is off,
	 * or when the specified category is not part of this ontology
	 */
	public void addSourceCategoryMapping(String edgeType, String category, boolean overwrite) throws IllegalArgumentException
	{
		if (!allSourceCategories.containsKey(category.toLowerCase()))
		{
			String errormsg = "The provided edge category ('" + category + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!overwrite && mapSourceTypeToCategory.containsKey(edgeType.toLowerCase()))
		{
			String errormsg = "The provided edge type is already mapped to a category!";
			throw new IllegalArgumentException(errormsg);
		}
		mapSourceTypeToCategory.put(edgeType.toLowerCase(), category.toLowerCase());
	}

	/**
	 * Remove all categories and type-categoryMappings
	 */
	protected void removeAllCategoriesAndMappings()
	{
		allSourceCategories = new HashMap<String, Boolean>();
		allSourceCategories.put(VOID_TYPE, true);

		allDiffCategories = new HashSet<String>();
		allDiffCategories.add(VOID_TYPE);

		mapSourceTypeToCategory = new HashMap<String, String>();
		addSourceCategoryMapping(VOID_TYPE, VOID_TYPE, false);
	}

	/**
	 * Remove all type-category mappings (keeping all defined categories).
	 */
	protected void removeAllSourceCategoryMappings()
	{
		mapSourceTypeToCategory = new HashMap<String, String>();
	}

}

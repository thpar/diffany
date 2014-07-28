package be.svlandeg.diffany.core.semantics;

import java.util.*;

import be.svlandeg.diffany.core.networks.*;
import be.svlandeg.diffany.core.visualstyle.EdgeDrawing;

/**
 * This class takes care of the semantic interpretation of different edge types and their corresponding categories in the 'source' networks,
 * i.e. the reference and condition-specific networks used as input.
 * It can define a differential edge by inspecting the two original edges (implemented by overriding classes).
 * Additionally, an overlapping edge can be generated by taking into account the edge type-to-category mapping.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class EdgeOntology
{

	private final String VOID_SYMMETRICAL_CAT;
	private final String VOID_DIRECTED_CAT;

	// edge types are stored in their canonical version (i.e. no punctuation, ascii-only)
	// when querying this map, always run the key through getCanonicalForm(edgeType) !
	protected Map<String, String> mapCanSourceTypeToCategory;

	protected Map<String, Boolean> allSourceCategories; // categories and their symmetry
	protected Set<String> posSourceCats;
	protected Set<String> negSourceCats;

	protected static Boolean default_symmetry = false;

	/**
	 * Create an empty ontology, with no edge categories or any mapping defined
	 */
	public EdgeOntology()
	{
		VOID_SYMMETRICAL_CAT = new EdgeGenerator().getVoidEdge(true).getType();
		VOID_DIRECTED_CAT = new EdgeGenerator().getVoidEdge(false).getType();
		
		removeAllCategoriesAndMappings();
		posSourceCats = new HashSet<String>();
		negSourceCats = new HashSet<String>();
	}
	
	/**
	 * Retrieve a void edge type
	 * 
	 * @param symmetrical whether or not the edge type should be symmetrical 
	 * @return the void edge type
	 */
	public String getVoidType(boolean symmetrical)
	{
		if (symmetrical)
		{
			return getCanonicalForm(VOID_SYMMETRICAL_CAT);
		}
		return getCanonicalForm(VOID_DIRECTED_CAT);
	}
	
	/**
	 * Retrieve a void edge category
	 * 
	 * @param symmetrical whether or not the edge category should be symmetrical 
	 * @return the void edge category
	 */
	public String getVoidCategory(boolean symmetrical)
	{
		if (symmetrical)
		{
			return VOID_SYMMETRICAL_CAT;
		}
		return VOID_DIRECTED_CAT;
	}

	/**
	 * Retrieve all the elements in this ontology that are root (i.e. do not have a parent)
	 * 
	 * @return all root elements in this ontology
	 */
	public abstract Set<String> retrieveAllSourceRootCats();

	/**
	 * Determine whether or not an edge type (child) is related to a category (partent).
	 * 
	 * @param childType the subclass edge type
	 * @param parentCat the superclass category
	 * @return whether or not the parent relationship holds, expressed by depth (-1 if unrelated, 0 if equal)
	 */
	public abstract int isSourceTypeChildOf(String childType, String parentCat);
	
	/**
	 * Determine whether or not two categories are related to eachother as child (sub) - parent (super)
	 * 
	 * @param childCat the subclass category
	 * @param parentCat the superclass category
	 * @return whether or not the parent relationship holds, expressed by depth (-1 if unrelated, 0 if equal)
	 */
	public abstract int isSourceCatChildOf(String childCat, String parentCat);
	
	/**
	 * Retrieve the parent category of a specific child category, or null if there is none.
	 * This method only goes one level up, so no grandparents etc. will be included.
	 * 
	 * @param childCat the subclass category
	 * @return the superclass category, or null if there is none
	 */
	public abstract String retrieveCatParent(String childCat);

	/**
	 * Return the common parent of two categories, or null if there is none
	 * 
	 * @param childCat1 the first child (sub) category
	 * @param childCat2 the second child (sub) category
	 * @return the common parent (super) category, or null if there is none such
	 */
	public abstract String commonSourceCatParent(String childCat1, String childCat2);
	
	/**
	 * For a set of categories, determine their most specific common parent/ancestor.
	 * Most specific is seen as a minimal maximum distance up to that ancestor across the whole categories set.
	 * 
	 * @param cats the original set of categories
	 * @return the most specific common parent, or null if there is none
	 */
	public abstract String retrieveFirstCommonParent(Collection<String> cats);

	/**
	 * Retrieve an {@link EdgeDrawing} object which knows how to define the visual styles in a differential network
	 * 
	 * @return the EdgeDrawing object for the differential network(s)
	 */
	public abstract EdgeDrawing getDifferentialEdgeDrawing();

	/**
	 * Retrieve an {@link EdgeDrawing} object which knows how to define the visual styles in an input (reference or condition-dependent) or overlapping network
	 * 
	 * @return the EdgeDrawing object for the input and overlapping network(s) (i.e. all but the differential networks)
	 */
	public abstract EdgeDrawing getSourceEdgeDrawing();

	/**
	 * Define a 'positive' category, like positive regulation
	 * 
	 * @param pos_cat the 'positive' category
	 * @throws IllegalArgumentException when the category is not defined in this ontology
	 */
	public void addPosSourceCat(String pos_cat) throws IllegalArgumentException
	{
		if (!isDefinedSourceCat(pos_cat))
		{
			String errormsg = "The positive category " + pos_cat + " is not defined in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		posSourceCats.add(pos_cat);
	}

	/**
	 * Define a 'negative' category, like negative regulation
	 * 
	 * @param neg_cat the 'negative' category
	 * @throws IllegalArgumentException when the category is not defined in this ontology
	 */
	public void addNegSourceCat(String neg_cat) throws IllegalArgumentException
	{
		if (!isDefinedSourceCat(neg_cat))
		{
			String errormsg = "The negative category is not defined in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		negSourceCats.add(neg_cat);
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
	 * Matching is done by comparing the canonical forms of the edge types (lowercase, alphanumerical only).
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return the semantic category of that edge or null if it is not mapped in this ontology
	 */
	public String getSourceCategory(String edgeType)
	{
		return mapCanSourceTypeToCategory.get(getCanonicalForm(edgeType));
	}

	/**
	 * Return the canonical form of an edge type. This form is obtained by transforming to lower case, removing all punctuation characters
	 * and retaining only alphanumerical characters. This enables straightforward matching of spelling variants.
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return the canonical form of the edgeType
	 */
	protected String getCanonicalForm(String edgeType)
	{
		String result = edgeType.toLowerCase();
		result = result.replaceAll("[^a-z0-9]+", "");
		return result;
	}

	/**
	 * Check whether a certain edge type is present in this ontology.
	 * Matching is done by comparing the canonical forms of the edge types (lowercase, alphanumerical only).
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return whether or not this edge type is present in this ontology
	 */
	public boolean isDefinedSourceType(String edgeType)
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
	 * Matching is done by comparing the canonical forms of the edge types (lowercase, alphanumerical only).
	 * 
	 * @param category the original category of the edge in a network
	 * @return whether or not this edge category is present in this ontology (will always be false when the category is null)
	 */
	public boolean isDefinedSourceCat(String category)
	{
		if (category == null)
		{
			return false;
		}
		return allSourceCategories.containsKey(category);
	}

	/**
	 * Return all source (input) categories present in this ontology.
	 * 
	 * @return the set of categories mapped in this ontology.
	 */
	public Set<String> getAllSourceCategories()
	{
		return allSourceCategories.keySet();
	}
	
	/**
	 * Return all positive source (input) categories present in this ontology.
	 * 
	 * @return the set of positive categories mapped in this ontology.
	 */
	public Set<String> getAllPosSourceCategories()
	{
		return posSourceCats;
	}
	
	/**
	 * Return all negative source (input) categories present in this ontology.
	 * 
	 * @return the set of negative categories mapped in this ontology.
	 */
	public Set<String> getAllNegSourceCategories()
	{
		return negSourceCats;
	}

	/**
	 * Retrieve the symmetry state of a source edge category in this ontology
	 * 
	 * @param edgeCat the source edge category
	 * @return the symmetry state of the source category
	 */
	public boolean isSymmetricalSourceCat(String edgeCat)
	{
		if (edgeCat == null || !allSourceCategories.containsKey(edgeCat))
		{
			String errormsg = "The source category should not be null or undefined!";
			throw new IllegalArgumentException(errormsg);
		}
		return allSourceCategories.get(edgeCat);
	}

	/**
	 * Retrieve the symmetry state of a source edgeType in this ontology
	 * 
	 * @param edgeType the source edgeType
	 * @return the symmetry state of the corresponding source category
	 */
	public boolean isSymmetricalSourceType(String edgeType)
	{
		if (edgeType == null || !isDefinedSourceType(edgeType))
		{
			return default_symmetry;
		}
		String cat = getSourceCategory(edgeType);
		return allSourceCategories.get(cat);
	}

	/**
	 * Add a source category
	 * 
	 * @param category the category that should be added to this ontology
	 * @param symmetric whether or not the category is symmetric
	 * @param overwrite determines whether or not this function may overwrite previous mappings of the same edge type
	 * 
	 * @throws IllegalArgumentException when the category is null, or already previously defined and overwrite is not switched on
	 */
	public void addSourceCategory(String category, boolean symmetric, boolean overwrite) throws IllegalArgumentException
	{
		if (category == null)
		{
			String errormsg = "The source category should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (! overwrite && allSourceCategories.containsKey(category))
		{
			String errormsg = "The source category " + category + " was already previously defined!";
			throw new IllegalArgumentException(errormsg);
		}
		allSourceCategories.put(category, symmetric);
	}

	/**
	 * Add a number of source categories (casing independent)
	 * 
	 * @param categories the categories that should be added to this ontology
	 * @param overwrite determines whether or not this function may overwrite previous mappings of the same edge type
	 */
	protected void addSourceCategories(Map<String, Boolean> categories, boolean overwrite)
	{
		for (String c : categories.keySet())
		{
			addSourceCategory(c, categories.get(c), overwrite);
		}
	}

	/**
	 * Create a new mapping from edge type to category. 
	 * Matching is done independent of upper/lower casing, considering only alphanumerical characters.
	 * 
	 * @param edgeType the original edge type - its canonical form should not have been defined in this ontology before
	 * @param category the category to be assigned to this edge type - this category should already exist in this ontology!
	 * @param overwrite determines whether or not this function may overwrite previous mappings of the same edge type
	 * 
	 * @throws IllegalArgumentException when the canonical form of the edge type was already mapped in this ontology and overwrite is off,
	 * or when the specified category is not part of this ontology
	 */
	public void addSourceCategoryMapping(String edgeType, String category, boolean overwrite) throws IllegalArgumentException
	{
		String canType = getCanonicalForm(edgeType);
		if (!allSourceCategories.containsKey(category))
		{
			String errormsg = "The provided edge category ('" + category + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!overwrite && mapCanSourceTypeToCategory.containsKey(canType))
		{
			String errormsg = "The provided edge type " + edgeType + "is already mapped to a category!";
			throw new IllegalArgumentException(errormsg);
		}
		mapCanSourceTypeToCategory.put(canType, category);
	}

	/**
	 * Remove all categories and type-categoryMappings
	 */
	protected void removeAllCategoriesAndMappings()
	{
		allSourceCategories = new HashMap<String, Boolean>();
		allSourceCategories.put(VOID_SYMMETRICAL_CAT, true);
		allSourceCategories.put(VOID_DIRECTED_CAT, false);

		mapCanSourceTypeToCategory = new HashMap<String, String>();
		addSourceCategoryMapping(VOID_SYMMETRICAL_CAT, VOID_SYMMETRICAL_CAT, false);
		addSourceCategoryMapping(VOID_DIRECTED_CAT, VOID_DIRECTED_CAT, false);
	}

}

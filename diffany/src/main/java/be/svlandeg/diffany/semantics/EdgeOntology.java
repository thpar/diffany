package be.svlandeg.diffany.semantics;

import java.awt.Paint;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This class takes care of the semantic interpretation of different edge types and their corresponding categories
 * in the 'source' networks, i.e. the reference and condition-specific networks used as input.
 * It can define a differential edge by inspecting the two original edges (implemented by overriding classes).
 * Additionally, an overlapping edge can be generated by taking into account the edge type-to-category mapping
 * and the direct super classes of categories.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class EdgeOntology
{

	public static final String VOID_TYPE = EdgeDefinition.getVoidEdge().getType().toLowerCase();

	private Map<String, String> sourceCatHierarchy;

	private Map<String, String> mapSourceTypeToCategory;
	private Set<String> allSourceCategories;

	private Set<String> allDiffCategories;

	/**
	 * Create an empty ontology, with no edge categories or any mapping defined
	 */
	public EdgeOntology()
	{
		removeAllCategoriesAndMappings();
		sourceCatHierarchy = new HashMap<String, String>();
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
	public EdgeDefinition getOverlapEdge(Set<EdgeDefinition> edges, double cutoff, boolean minOperator) throws IllegalArgumentException
	{
		if (edges == null || edges.isEmpty())
		{
			String errormsg = "The set of edges should not be null or empty!";
			throw new IllegalArgumentException(errormsg);
		}
		EdgeDefinition overlap_edge = new EdgeDefinition();
		int countEdges = edges.size();
		
		//////////// DETERMINE NEGATION AND SYMMETRY ////////////////
		int countNegated = 0;
		int countSymmetrical = 0;

		double minWeight = Double.MAX_VALUE;
		double maxWeight = Double.MIN_VALUE;

		for (EdgeDefinition e : edges)
		{
			if (e.isNegated())
				countNegated++;

			if (e.isSymmetrical())
				countSymmetrical++;

			double weight = e.getWeight();
			if (weight < minWeight)
				minWeight = weight;
			if (weight > maxWeight)
				maxWeight = weight;
		}

		 // some are negated, some are not -> no overlap
		if (countNegated != 0 && countNegated != countEdges)
		{
			return EdgeDefinition.getVoidEdge();
		}
		boolean overlapNegated = countNegated == countEdges;
		overlap_edge.makeNegated(overlapNegated);

		overlap_edge.makeSymmetrical(countSymmetrical == countEdges);
		
		//////////// DETERMINE WEIGHT ////////////////

		// the overlapping weight is the minimum between the two, or the maximum
		// if specified as such
		double overlapWeight = minWeight;
		if (!minOperator)
		{
			overlapWeight = maxWeight;
		}
		if (overlapWeight <= cutoff)
		{
			return EdgeDefinition.getVoidEdge();
		}
		overlap_edge.setWeight(overlapWeight);
		
		//////////// DEFINE TYPE BY INSPECTING CHILDREN AND PARENTS ////////////////

		Map<String, Integer> allCommonParents = retrieveFirstCommonParents(edges, false);
		if (allCommonParents.isEmpty()) 
		{
			// no category covers all of the edges
			return EdgeDefinition.getVoidEdge();
		}
		String commonParent = allCommonParents.keySet().iterator().next();
		int minParentDepth = allCommonParents.get(commonParent);
		
		
		if (minParentDepth == 0 || !overlapNegated)	//  the shared edge is the (first) common super class 
		{
			// if there are multiple common parents, equally close (minimal depth), we take the first one at random
			overlap_edge.setType(commonParent);
			return overlap_edge;
		} 
		
		Map<String, Integer> allCommonChildren = retrieveFirstCommonChildren(edges);
		
		if (! allCommonChildren.isEmpty()) 	
		{
			// the shared edge is the negation of the (first) common subclass, if there is one such
			String commonChild = allCommonChildren.keySet().iterator().next();
			overlap_edge.setType(commonChild);
			return overlap_edge;
		}

		return EdgeDefinition.getVoidEdge();
	}
	
	/**
	 * For a set of EdgeDefinition objects, determine their most general common children.
	 * Most general is seen as a minimal maximum distance down to that (grand)child across the whole edge set.
	 * @param edges the original set of edges
	 * @return a map of most general common children and their (equal) maximal distance to the original edges.
	 */
	protected Map<String, Integer> retrieveFirstCommonChildren(Set<EdgeDefinition> edges)
	{
		int countEdges = edges.size();
		Map<String, Integer> allCommonChildren = new HashMap<String, Integer>();
		
		if (edges.size() == 1)
		{
			String cat = getSourceCategory(edges.iterator().next().getType());
			allCommonChildren.put(cat, 0);
			return allCommonChildren;
		}
		
		// go down the ontology tree and fetch all  children
		Map<String, Integer> childrenCatsByCount = new HashMap<String, Integer>();
				
		// count the depth (down) in the ontology tree
		Map<String, Integer> childByDepth = new HashMap<String, Integer>();

		for (EdgeDefinition e : edges)
		{
			int depth = 0;
			String cat = getSourceCategory(e.getType());
			// each cat is its own child of depth 0
			addOne(childrenCatsByCount, cat);
			recordMaxDepth(childByDepth, cat, depth);
					
			// record all children down in the hierarchy
			Set<String> childrenCats = retrieveChildren(cat);
			while (childrenCats.size() > 0)
			{
				Set<String> newChildrenCats = new HashSet<String>();
				depth++;
				for (String child : childrenCats)
				{
					addOne(childrenCatsByCount, child);
					recordMaxDepth(childByDepth, child, depth);
					newChildrenCats.addAll(retrieveChildren(child));
				}
				childrenCats = newChildrenCats;
			}
		}
		
		// find the closest common child (if any)
		for (String cat : childrenCatsByCount.keySet())
		{
			if (childrenCatsByCount.get(cat) == countEdges)
				allCommonChildren.put(cat, childByDepth.get(cat));
		}
		int minChildDepth = Integer.MAX_VALUE;
		for (String foundCat : allCommonChildren.keySet())
		{
			int foundDepth = childByDepth.get(foundCat);
			minChildDepth = Math.min(minChildDepth, foundDepth);
		}
		for (String cat : childByDepth.keySet())
		{
			int thisDepth = childByDepth.get(cat);
			if (thisDepth != minChildDepth)
			{
				// don't keep more specific categories if there is a more general common child
				allCommonChildren.remove(cat);	
			}
		}
		
		return allCommonChildren;
	}
	
	/**
	 * For a set of EdgeDefinition objects, determine their most specific common parents/ancestors.
	 * Most specific is seen as a minimal maximum distance up to that ancestor across the whole edge set.
	 * @param edges the original set of edges
	 * @return a map of most specific common parents and their (equal) maximal distance to the original edges.
	 */
	protected Map<String, Integer> retrieveFirstCommonParents(Set<EdgeDefinition> edges, boolean excludeEmpty)
	{
		int countEdges = edges.size();
		Map<String, Integer> allCommonParents = new HashMap<String, Integer>();
		
		if (countEdges == 1)
		{
			String cat = getSourceCategory(edges.iterator().next().getType());
			allCommonParents.put(cat, 0);
			return allCommonParents;
		}
		
		// go up the ontology tree and fetch all parents
		Map<String, Integer> parentCatsByCount = new HashMap<String, Integer>();
				
		// count the depth (up) in the ontology tree
		Map<String, Integer> parentByDepth = new HashMap<String, Integer>();
		
		int countEmpty = 0;

		for (EdgeDefinition e : edges)
		{
			int depth = 0;
			String cat = getSourceCategory(e.getType());
			if (cat.equals(VOID_TYPE))
			{
				countEmpty++;
			}
			else
			{
				// each cat is its own parent of depth 0
				addOne(parentCatsByCount, cat);	
				recordMaxDepth(parentByDepth, cat, depth);
					
				// record all parents up in the hierarchy
				String parentCat = retrieveParent(cat);
				while (parentCat != null)
				{
					depth++;
					addOne(parentCatsByCount, parentCat);
					recordMaxDepth(parentByDepth, parentCat, depth);
					parentCat = retrieveParent(parentCat);
				}
			}
		}
		if (countEmpty == countEdges)
		{
			allCommonParents.put(VOID_TYPE, 0);
			return allCommonParents;
		}
		if (excludeEmpty)
			countEdges = countEdges - countEmpty;
		
		for (String cat : parentCatsByCount.keySet())
		{
			if (parentCatsByCount.get(cat) == countEdges)
				allCommonParents.put(cat, parentByDepth.get(cat));
		}
		if (allCommonParents.keySet().size() < 1) 
		{
			// no category covers all of the edges
			return new HashMap<String, Integer>();
		}
		
		int minParentDepth = Integer.MAX_VALUE;
		for (String foundCat : allCommonParents.keySet())
		{
			int foundDepth = parentByDepth.get(foundCat);
			minParentDepth = Math.min(minParentDepth, foundDepth);
		}
		for (String cat : parentByDepth.keySet())
		{
			int thisDepth = parentByDepth.get(cat);
			if (thisDepth != minParentDepth)
			{
				// don't keep more general categories if there is a more specific common parent
				allCommonParents.remove(cat);	
			}
		}
		return allCommonParents;
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
		return allSourceCategories.contains(category);
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
	 * Define the visual style of an edge in a differential network, by edge category.
	 * 
	 * @param category the category of the edge interaction
	 * @return a Paint object which specifies how the edge should be drawn
	 */
	public abstract Paint getDifferentialEdgeStyle(String category);
	
	/**
	 * Define the visual style of an edge in a 'normal' network (reference, condition-dependent or overlap).
	 * 
	 * @param edgeType the type of the edge interaction
	 * @return a Paint object which specifies how the edge should be drawn
	 */
	public abstract Paint getSourceEdgeStyle(String edgeType);


	/**
	 * Return all source (input) categories present in this ontology.
	 * @return the set of categories mapped in this ontology.
	 */
	public Set<String> getAllSourceCategories()
	{
		return allSourceCategories;
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
			String errormsg = "The category should not be null!";
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
			allDiffCategories.add(c.toLowerCase());
		}
	}

	/**
	 * Add a source category (casing independent)
	 * @param category the category that should be added to this ontology
	 * @throws IllegalArgumentException when the category is null
	 */
	protected void addSourceCategory(String category) throws IllegalArgumentException
	{
		if (category == null)
		{
			String errormsg = "The category should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		allSourceCategories.add(category.toLowerCase());
	}

	/**
	 * Add a number of source categories (casing independent)
	 * @param categories the categories that should be added to this ontology
	 */
	protected void addSourceCategories(Set<String> categories)
	{
		for (String c : categories)
		{
			allSourceCategories.add(c.toLowerCase());
		}
	}

	/**
	 * Define a child category and its parent category. The child should not have received a parent before.
	 * @param childCat the child category (subclass)
	 * @param parentCat the parent category (superclass)
	 * @throws IllegalArgumentException when the childCat was already previously
	 * attached to a parent or if either of the two categories are not defined in this ontology
	 */
	protected void putSourceParent(String childCat, String parentCat) throws IllegalArgumentException
	{
		if (sourceCatHierarchy.containsKey(childCat))
		{
			String errormsg = "The provided child category ('" + childCat + "') already has a parent category!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!allSourceCategories.contains(childCat.toLowerCase()))
		{
			String errormsg = "The provided child category ('" + childCat + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!allSourceCategories.contains(parentCat.toLowerCase()))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		sourceCatHierarchy.put(childCat, parentCat);
	}
	
	
	/**
	 * Retrieve the parent category of a specific child category, or null if there is none.
	 * @param childCat the subclass category
	 * @return the superclass category, or null if there is none
	 */
	protected String retrieveParent(String childCat)
	{
		return sourceCatHierarchy.get(childCat);
	}
	
	/**
	 * Retrieve the set of child categories of a specific parent category, or an empty set if there are none.
	 * @param parentCat the superclass category
	 * @return the set of subclass categories, or an empty set if there are none.
	 */
	protected Set<String> retrieveChildren(String parentCat)
	{
		Set<String> children = new HashSet<String>();
		for (String childCat : sourceCatHierarchy.keySet())
		{
			if (sourceCatHierarchy.get(childCat).equals(parentCat))
				children.add(childCat);
		}
		return children;
	}

	/**
	 * Determine whether or not two categories are related to eachother as child (sub) - parent (super)
	 * @param childCat the subclass category
	 * @param parentCat the superclass category
	 * @return whether or not the parent relationship holds, expressed by depth (-1 if unrelated, 0 if equal)
	 */
	protected int isSourceChildOf(String childCat, String parentCat)
	{
		int depth = 0;
		if (childCat.equals(parentCat))
		{
			return 0;
		}
		String parent = retrieveParent(childCat);
		while (parent != null)
		{
			depth++;
			if (parent.equals(parentCat))
			{
				return depth;
			}
			parent = retrieveParent(parent);
		}
		return -1;
	}

	/**
	 * Return the common parent of two categories, or null if there is none
	 * @param childCat1 the first child (sub) category
	 * @param childCat2 the second child (sub) category
	 * @return the common parent (super) category, or null if there is none such
	 */
	protected String commonSourceParent(String childCat1, String childCat2)
	{
		if (!sourceCatHierarchy.containsKey(childCat1))
		{
			return null;
		}
		if (!sourceCatHierarchy.containsKey(childCat2))
		{
			return null;
		}
		if (sourceCatHierarchy.get(childCat1).equals(sourceCatHierarchy.get(childCat2)))
		{
			return sourceCatHierarchy.get(childCat1);
		}
		return null;
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
		if (!allSourceCategories.contains(category.toLowerCase()))
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
		allSourceCategories = new HashSet<String>();
		allSourceCategories.add(VOID_TYPE);

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

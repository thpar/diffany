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
	
	private Map<String, String> catHierarchy;

	private Map<String, String> mapEdgeToCategory;
	private Set<String> allCategories;

	/**
	 * Create an empty ontology, with no edge categories or any mapping defined
	 */
	public EdgeOntology()
	{
		removeAllCategoriesAndMappings();
		catHierarchy = new HashMap<String, String>();
	}

	/**
	 * Method that defines the differential edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns EdgeDefinition.getVoidEdge() when the edge should be deleted (i.e. not present in differential network).
	 * 
	 * @param refEdge the edge definition in the reference network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param conEdge the edge definition in the condition-specific network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential network
	 * 
	 * @return the edge definition in the differential network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public abstract EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, EdgeDefinition conEdge, double cutoff) throws IllegalArgumentException;

	/**
	 * Method that defines the overlapping edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns EdgeDefinition.getVoidEdge() when the edge should be deleted (i.e. not present in the shared network).
	 * 
	 * @param refEdge the edge definition in the reference network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param conEdge the edge definition in the condition-specific network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param cutoff the minimal value of a resulting edge for it to be included in the shared network
	 * 
	 * @return the edge definition in the shared network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public EdgeDefinition getSharedEdge(EdgeDefinition refEdge, EdgeDefinition conEdge, double cutoff) throws IllegalArgumentException
	{
		EdgeDefinition shared_edge = new EdgeDefinition();

		String refCat = getCategory(refEdge.getType());
		String conCat = getCategory(conEdge.getType());
		
		boolean refNeg = refEdge.isNegated();
		boolean conNeg = conEdge.isNegated();

		if (refNeg != conNeg)
		{
			return EdgeDefinition.getVoidEdge();
		}
		
		shared_edge.makeNegated(refNeg);
			
		// the shared weight is the minimum between the two
		double sharedWeight = Math.min(refEdge.getWeight(), conEdge.getWeight());
		if (sharedWeight <= cutoff)
		{
			return EdgeDefinition.getVoidEdge();
		}
		shared_edge.setWeight(sharedWeight);

		// the shared edge is only symmetrical if both original edges are
		boolean refSymm = refEdge.isSymmetrical();
		boolean conSymm = conEdge.isSymmetrical();
		boolean sharedSymm = refSymm && conSymm;
		shared_edge.makeSymmetrical(sharedSymm);

		if (refCat.equals(conCat))
		{
			shared_edge.setType(refEdge.getType());
			return shared_edge;
		}
		if (isChildOf(refCat,conCat))	// conCat is superclass
		{
			if (conNeg)	// there is negation -> shared edge is negated of type subclass
			{
				shared_edge.setType(refEdge.getType());
				return shared_edge;
			}
			else	// there is no negation -> shared edge is of type superclass
			{
				shared_edge.setType(conEdge.getType());
				return shared_edge;
			}
		}
		if (isChildOf(conCat,refCat))	// refCat is superclass
		{
			if (refNeg)	// there is negation -> shared edge is negated of type subclass
			{
				shared_edge.setType(conEdge.getType());
				return shared_edge;
			}
			else	// there is no negation -> shared edge is of type superclass
			{
				shared_edge.setType(refEdge.getType());
				return shared_edge;
			}
		}
		
		return EdgeDefinition.getVoidEdge();
		
	}
	
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
	 * Check whether a certain edge type is present in this ontology. Matching is done independent of upper/lower casing.
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return whether or not this edge type is present in this ontology
	 */
	protected boolean isDefined(String edgeType)
	{
		if (edgeType == null)
		{
			return false;
		}
		return mapEdgeToCategory.containsKey(edgeType.toLowerCase());
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
	 * Define a child category and its parent category. The child should not have received a parent before.
	 * @param childCat the child category (subclass)
	 * @param parentCat the parent category (superclass)
	 * @throws IllegalArgumentException when the childCat was already previously attached to a parent or if either of the two categories are not defined in this ontology
	 */
	protected void putParent(String childCat, String parentCat) throws IllegalArgumentException
	{
		if (catHierarchy.containsKey(childCat))
		{
			String errormsg = "The provided child category ('" + childCat + "') already has a parent category!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!allCategories.contains(childCat.toLowerCase()))
		{
			String errormsg = "The provided child category ('" + childCat + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!allCategories.contains(parentCat.toLowerCase()))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		catHierarchy.put(childCat, parentCat);
	}
	
	/**
	 * Determine whether or not two categories are related to eachother as child (sub) - parent (super)
	 * @param childCat the subclass category
	 * @param parentCat the superclass category
	 * @return whether or not the parent relationship holds
	 */
	protected boolean isChildOf(String childCat, String parentCat)
	{
		if (! catHierarchy.containsKey(childCat))
		{
			return false;
		}
		return catHierarchy.get(childCat).equals(parentCat);
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

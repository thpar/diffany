package be.svlandeg.diffany.core.semantics;

import java.util.*;

import be.svlandeg.diffany.core.networks.EdgeDefinition;

/**
 * This edge ontology is implemented as a tree structure.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class TreeEdgeOntology extends EdgeOntology
{

	// map children to parents
	protected Map<String, String> sourceCatHierarchy;

	// prefixes for easy classification of types
	private String negPrefix_symm;
	private String negPrefix_dir;
	private String posPrefix_symm;
	private String posPrefix_dir;
	private String unspecifiedPrefix;

	/**
	 * Create a new ontology, defining the set of categories. a
	 * After the constructor is called, default edge-category mappings should be inserted!
	 * 
	 * @param posPrefix_symm the prefix to be used when a symmetrical interaction increases
	 * @param posPrefix_dir the prefix to be used when a directed interaction increases
	 * @param negPrefix_symm the prefix to be used when a symmetrical interaction decreases
	 * @param negPrefix_dir the prefix to be used when a directed interaction decreases
	 * @param unspecifiedPrefix the prefix to be used when an interaction is unspecific
	 * @param sourceCats all the source categories defined in this ontology
	 */
	public TreeEdgeOntology(String posPrefix_symm, String posPrefix_dir, String negPrefix_symm, String negPrefix_dir, String unspecifiedPrefix, Map<String, Boolean> sourceCats)
	{
		super();
		sourceCatHierarchy = new HashMap<String, String>();
		this.negPrefix_symm = negPrefix_symm;
		this.negPrefix_dir = negPrefix_dir;
		this.posPrefix_symm = posPrefix_symm;
		this.posPrefix_dir = posPrefix_dir;
		this.unspecifiedPrefix = unspecifiedPrefix;

		addSourceCategories(sourceCats, false);
	}

	////////////// TREE ALGORITHMS //////////////////////////////////

	@Override
	public Set<String> retrieveAllSourceRootCats(boolean ignoreGenerics)
	{
		Set<String> allRoots = new HashSet<String>();
		for (String cat : allSourceCategories.keySet())
		{
			String parent = sourceCatHierarchy.get(cat);
			if (parent == null) 
			{
				allRoots.add(cat);
			}
			else if (ignoreGenerics)
			{
				if (parent.equals(GENERIC_SYMMETRICAL_CAT) || parent.equals(GENERIC_DIRECTED_CAT))
				{
					allRoots.add(cat);
				}
			}
		}
		if (ignoreGenerics)
		{
			allRoots.remove(GENERIC_SYMMETRICAL_CAT);
			allRoots.remove(GENERIC_DIRECTED_CAT);
		}
		return allRoots;
	}

	/**
	 * Define a child category and its parent category. The child should not have received a parent before.
	 * 
	 * @param childCat the child category (subclass)
	 * @param parentCat the parent category (superclass)
	 * @throws IllegalArgumentException when the childCat was already previously
	 * attached to a parent or if either of the two categories are not defined in this ontology
	 */
	protected void putSourceCatParent(String childCat, String parentCat) throws IllegalArgumentException
	{
		if (sourceCatHierarchy.containsKey(childCat))
		{
			String errormsg = "The provided child category ('" + childCat + "') already has a parent category!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!allSourceCategories.containsKey(parentCat))
		{
			String errormsg = "The provided child category ('" + childCat + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!allSourceCategories.containsKey(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		sourceCatHierarchy.put(childCat, parentCat);
	}

	@Override
	public String retrieveCatParent(String childCat)
	{
		return sourceCatHierarchy.get(childCat);
	}

	/**
	 * Retrieve the set of child categories of a specific parent category, or an empty set if there are none and this category is thus a 'leaf'.
	 * This method only goes one level deep, so no grandchildren etc. will be included.
	 * 
	 * @param parentCat the superclass category
	 * @return the set of subclass categories, or an empty set if there are none
	 */
	public Set<String> retrieveCatChildren(String parentCat)
	{
		Set<String> children = new HashSet<String>();
		for (String childCat : sourceCatHierarchy.keySet())
		{
			if (sourceCatHierarchy.get(childCat).equals(parentCat))
			{
				children.add(childCat);
			}
		}
		return children;
	}

	@Override
	public int isSourceTypeChildOf(String childType, String grandparentCat)
	{
		if (childType == null || grandparentCat == null)
		{
			return -1;
		}
		String cat = getSourceCategory(childType);
		return isSourceCatChildOf(cat, grandparentCat);
	}

	@Override
	public int isSourceCatChildOf(String childCat, String grandparentCat)
	{
		if (childCat == null || grandparentCat == null)
		{
			return -1;
		}
		int depth = 0;

		if (childCat.equals(grandparentCat))
		{
			return 0;
		}

		String parent = retrieveCatParent(childCat);

		while (parent != null)
		{
			depth++;
			if (parent.equals(grandparentCat))
			{
				return depth;
			}
			parent = retrieveCatParent(parent);
		}
		return -1;
	}

	@Override
	public String commonSourceCatParent(String childCat1, String childCat2)
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
	 * For a set of EdgeDefinition objects, determine all their common children.
	 * 
	 * @param edges the original set of edges
	 * @return a map of all common children and their maximal distance to the original edges.
	 */
	protected Map<String, Integer> retrieveAllCommonChildren(Collection<EdgeDefinition> edges)
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
			Set<String> childrenCats = retrieveCatChildren(cat);
			while (childrenCats.size() > 0)
			{
				Set<String> newChildrenCats = new HashSet<String>();
				depth++;
				for (String child : childrenCats)
				{
					addOne(childrenCatsByCount, child);
					recordMaxDepth(childByDepth, child, depth);
					newChildrenCats.addAll(retrieveCatChildren(child));
				}
				childrenCats = newChildrenCats;
			}
		}

		// find the closest common child (if any)
		for (String cat : childrenCatsByCount.keySet())
		{
			if (childrenCatsByCount.get(cat) == countEdges)
			{
				allCommonChildren.put(cat, childByDepth.get(cat));
			}
		}

		return allCommonChildren;
	}

	/**
	 * For a set of EdgeDefinition objects, determine their most general common child.
	 * Most general is seen as a minimal maximum distance down to that (grand)child across the whole edge set.
	 * 
	 * @param edges the original set of edges
	 * @return the most general common child.
	 */
	public String retrieveFirstCommonChild(Collection<EdgeDefinition> edges)
	{
		Map<String, Integer> allCommonChildren = retrieveAllCommonChildren(edges);
		Set<String> allFirstChildren = new HashSet<String>();

		int minChildDepth = Integer.MAX_VALUE;
		for (String foundCat : allCommonChildren.keySet())
		{
			int foundDepth = allCommonChildren.get(foundCat);
			minChildDepth = Math.min(minChildDepth, foundDepth);
		}
		for (String cat : allCommonChildren.keySet())
		{
			int thisDepth = allCommonChildren.get(cat);
			if (thisDepth == minChildDepth)
			{
				// don't keep more specific categories if there is a more general common child
				allFirstChildren.add(cat);
			}
		}
		if (allFirstChildren.isEmpty())
		{
			return null;
		}
		if (allFirstChildren.size() > 1)
		{
			String errormsg = "Found more than 1 first common child!";
			throw new RuntimeException(errormsg);
		}
		return allFirstChildren.iterator().next();
	}

	/**
	 * For a set of edge types, determine all their common parents/ancestors.
	 * 
	 * @param cats the original set of edge types
	 * @return a map of all common parents and their maximal distance to the original edges
	 */
	protected Map<String, Integer> retrieveAllCommonParents(Collection<String> cats)
	{
		int countEdges = cats.size();
		Map<String, Integer> allCommonParents = new HashMap<String, Integer>();

		if (countEdges == 1)
		{
			String cat = cats.iterator().next();
			allCommonParents.put(cat, 0);
			return allCommonParents;
		}

		// go up the ontology tree and fetch all parents
		Map<String, Integer> parentCatsByCount = new HashMap<String, Integer>();

		// count the depth (up) in the ontology tree
		Map<String, Integer> parentByDepth = new HashMap<String, Integer>();

		for (String cat : cats)
		{
			int depth = 0;
			// each cat is its own parent of depth 0
			addOne(parentCatsByCount, cat);
			recordMaxDepth(parentByDepth, cat, depth);

			// record all parents up in the hierarchy
			String parentCat = retrieveCatParent(cat);
			while (parentCat != null)
			{
				depth++;
				addOne(parentCatsByCount, parentCat);
				recordMaxDepth(parentByDepth, parentCat, depth);
				parentCat = retrieveCatParent(parentCat);
			}

		}

		for (String cat : parentCatsByCount.keySet())
		{
			if (parentCatsByCount.get(cat) == countEdges)
			{
				allCommonParents.put(cat, parentByDepth.get(cat));
			}
		}
		if (allCommonParents.keySet().size() < 1)
		{
			// no category covers all of the edges
			return new HashMap<String, Integer>();
		}
		return allCommonParents;
	}

	@Override
	public String retrieveFirstCommonParent(Collection<String> cats)
	{
		Map<String, Integer> allCommonParents = retrieveAllCommonParents(cats);
		Set<String> allFirstParents = new HashSet<String>();

		int minParentDepth = Integer.MAX_VALUE;
		for (String foundCat : allCommonParents.keySet())
		{
			int foundDepth = allCommonParents.get(foundCat);
			minParentDepth = Math.min(minParentDepth, foundDepth);
		}
		for (String cat : allCommonParents.keySet())
		{
			int thisDepth = allCommonParents.get(cat);
			if (thisDepth == minParentDepth)
			{
				// don't keep more general categories if there is a more specific common parent
				allFirstParents.add(cat);
			}
		}
		if (allFirstParents.isEmpty())
		{
			return null;
		}
		if (allFirstParents.size() > 1)
		{
			String errormsg = "Found more than 1 first common parent!";
			throw new RuntimeException(errormsg);
		}
		return allFirstParents.iterator().next();
	}

	////////////// DIFFERENTIAL EDGE STATE //////////////////////////////////

	/**
	 * Retrieve whether or not a certain differential category can be seen as a positive, directed edge
	 * 
	 * @param category the category of the edge interaction
	 * @return whether or not a certain category can be seen as a positive, directed edge
	 */
	public boolean isPosDirected(String category)
	{
		return category.startsWith(posPrefix_dir);

	}

	/**
	 * Retrieve whether or not a certain differential category can be seen as a negative, directed edge
	 * 
	 * @param category the category of the edge interaction
	 * @return whether or not a certain differential category can be seen as a negative, directed edge
	 */
	public boolean isNegDirected(String category)
	{
		return category.startsWith(negPrefix_dir);

	}

	/**
	 * Retrieve whether or not a certain category can be seen as a negative, symmetrical edge
	 * 
	 * @param category the category of the edge interaction
	 * @return whether or not a certain category can be seen as a negative, symmetrical edge
	 */
	public boolean isNegSymm(String category)
	{
		return category.startsWith(negPrefix_symm);

	}

	/**
	 * Retrieve whether or not a certain category can be seen as a positive, symmetrical edge
	 * 
	 * @param category the category of the edge interaction
	 * @return whether or not a certain category can be seen as a positive, symmetrical edge
	 */
	public boolean isPosSymm(String category)
	{
		return category.startsWith(posPrefix_symm);

	}

	/**
	 * Retrieve the prefix used for negative symmetrical edges
	 * 
	 * @return the prefix used for negative symmetrical edges
	 */
	public String getNegPrefix_symm()
	{
		return negPrefix_symm;
	}

	/**
	 * Retrieve the prefix used for negative directed edges
	 * 
	 * @return the prefix used for negative directed edges
	 */
	public String getNegPrefix_dir()
	{
		return negPrefix_dir;
	}

	/**
	 * Retrieve the prefix used for positive symmetrical edges
	 * 
	 * @return the prefix used for positive symmetrical edges
	 */
	public String getPosPrefix_symm()
	{
		return posPrefix_symm;
	}

	/**
	 * Retrieve the prefix used for positive directed edges
	 * 
	 * @return the prefix used for positive directed edges
	 */
	public String getPosPrefix_dir()
	{
		return posPrefix_dir;
	}

	/**
	 * Retrieve the prefix used for unspecified edges
	 * 
	 * @return the prefix used for unspecified edges
	 */
	public String getUnspecifiedPrefix()
	{
		return unspecifiedPrefix;
	}

}

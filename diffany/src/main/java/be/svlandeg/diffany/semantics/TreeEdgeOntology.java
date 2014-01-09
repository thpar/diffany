package be.svlandeg.diffany.semantics;

import java.awt.Color;
import java.util.*;

import be.svlandeg.diffany.concepts.EdgeDefinition;
import be.svlandeg.diffany.concepts.VisualEdgeStyle.ArrowHead;

/**
 * This edge ontology is implemented as a tree structure.
 * 
 * @author Sofie Van Landeghem
 */
public class TreeEdgeOntology extends EdgeOntology
{
	
	protected Map<String, String> sourceCatHierarchy;

	private String negPrefix_symm;
	private String negPrefix_dir;
	private String posPrefix_symm;
	private String posPrefix_dir;
	
	protected static Color neg_diff_paint = Color.RED;
	protected static Color pos_diff_paint = Color.GREEN;
	protected static Color default_diff_paint = Color.GRAY;
	
	protected static ArrowHead neg_diff_ah = ArrowHead.ARROW;
	protected static ArrowHead pos_diff_ah = ArrowHead.ARROW;
	protected static ArrowHead default_diff_ah = ArrowHead.ARROW;
	
	protected static Color neutral_source_paint = Color.LIGHT_GRAY;
	protected static ArrowHead neutral_source_ah = ArrowHead.ARROW;
	
	protected static ArrowHead symm_ah = ArrowHead.NONE;
	
	protected Map<String, Color> parentSourceCatToColor;
	protected Map<String, ArrowHead> parentSourceCatToArrowHead;
	

	/**
	 * Create a new ontology, defining the set of categories. and inserting
	 * After the constructor is called, default edge-category mappings should be inserted!
	 * @param posPrefix_symm the prefix to be used when a symmetrical interaction increases
	 * @param posPrefix_dir the prefix to be used when a directed interaction increases
	 * @param negPrefix_symm the prefix to be used when a symmetrical interaction decreases
	 * @param negPrefix_dir the prefix to be used when a directed interaction decreases
	 * @param sourceCats all the source categories defined in this ontology
	 */
	public TreeEdgeOntology(String posPrefix_symm, String posPrefix_dir, String negPrefix_symm, String negPrefix_dir, Map<String, Boolean> sourceCats)
	{
		super();
		sourceCatHierarchy = new HashMap<String, String>();
		this.negPrefix_symm = negPrefix_symm;
		this.negPrefix_dir = negPrefix_dir;
		this.posPrefix_symm = posPrefix_symm;
		this.posPrefix_dir = posPrefix_dir;
		parentSourceCatToColor = new HashMap<String, Color>();
		parentSourceCatToArrowHead = new HashMap<String, ArrowHead>();
		addSourceCategories(sourceCats);
	}
	
	
	////////////// TREE ALGORITHMS //////////////////////////////////
	
	
	@Override
	public Set<String> retrieveAllSourceRootCats()
	{
		Set<String> allRoots = new HashSet<String>();
		for (String cat : allSourceCategories.keySet())
		{
			if (! sourceCatHierarchy.keySet().contains(cat))	// does not have a parent
			{
				allRoots.add(cat);
			}
		}
		return allRoots;
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
		if (!allSourceCategories.containsKey(childCat.toLowerCase()))
		{
			String errormsg = "The provided child category ('" + childCat + "') does not exist in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		if (!allSourceCategories.containsKey(parentCat.toLowerCase()))
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
	
	@Override
	public int isSourceChildOf(String childCat, String parentCat)
	{
		if (childCat == null || parentCat == null)
		{
			return -1;
		}
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

	@Override
	public String commonSourceParent(String childCat1, String childCat2)
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
		
		return allCommonChildren;
	}
	
	/**
	 * For a set of EdgeDefinition objects, determine their most general common child.
	 * Most general is seen as a minimal maximum distance down to that (grand)child across the whole edge set.
	 * @param edges the original set of edges
	 * @return the most general common child.
	 */
	protected String retrieveFirstCommonChild(Collection<EdgeDefinition> edges)
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
	 * For a set of EdgeDefinition objects, determine all their common parents/ancestors.
	 * 
	 * @param edges the original set of edges
	 * @param excludeEmpty whether or not to exclude empty edges when looking for a common ancestor
	 * @return a map of all common parents and their maximal distance to the original edges
	 */
	protected Map<String, Integer> retrieveAllCommonParents(Collection<EdgeDefinition> edges, boolean excludeEmpty)
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
		return allCommonParents;
	}
	
	/**
	 * For a set of EdgeDefinition objects, determine their most specific common parent/ancestor.
	 * Most specific is seen as a minimal maximum distance up to that ancestor across the whole edge set.
	 * 
	 * @param edges the original set of edges
	 * @param excludeEmpty whether or not to exclude empty edges when looking for a common ancestor
	 * @return the most specific common parent, or null if there is none
	 */
	protected String retrieveFirstCommonParent(Collection<EdgeDefinition> edges, boolean excludeEmpty)
	{
		Map<String, Integer> allCommonParents = retrieveAllCommonParents(edges, excludeEmpty);
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
	
	//////////////// ACTUAL EDGE ONTOLOGY METHODS ////////////////
	
	@Override
	public EdgeDefinition getOverlapEdge(Collection<EdgeDefinition> edges, double cutoff, boolean minOperator) throws IllegalArgumentException
	{
		if (edges == null || edges.isEmpty())
		{
			String errormsg = "The set of edges should not be null or empty!";
			throw new IllegalArgumentException(errormsg);
		}
		EdgeDefinition overlap_edge = new EdgeDefinition();
		int countEdges = edges.size();
		
		// 1. DETERMINE NEGATION AND SYMMETRY //
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
		
		// 2. DETERMINE WEIGHT //

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
		
		// 3. DEFINE TYPE BY INSPECTING CHILDREN AND PARENTS //

		String firstCommonParent = retrieveFirstCommonParent(edges, false);
		if (firstCommonParent == null) 
		{
			// no category covers all of the edges
			return EdgeDefinition.getVoidEdge();
		}
		
		if (!overlapNegated && firstCommonParent != null)	//  the shared edge is the (first) common super class 
		{
			overlap_edge.setType(firstCommonParent);
			return overlap_edge;
		} 
		
		String firstCommonChild = retrieveFirstCommonChild(edges);
		if (firstCommonParent == null) 
		{
			// no category covers all of the edges
			return EdgeDefinition.getVoidEdge();
		}
		
		// the shared edge is the negation of the (first) common subclass, if there is one such
		overlap_edge.setType(firstCommonChild);
		return overlap_edge;
	}

	
	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, Collection<EdgeDefinition> conEdges, double cutoff) throws IllegalArgumentException
	{
		EdgeDefinition diff_edge = new EdgeDefinition();
		Set<EdgeDefinition> conEdges2 = new HashSet<EdgeDefinition>();
		Set<EdgeDefinition> allEdges = new HashSet<EdgeDefinition>();
		
		//////////// DEAL WITH SYMMETRY AND NEGATION ////////////////
		boolean conSymm = true;
		for (EdgeDefinition e : conEdges)
		{
			// negated edges are set to void
			if (e.isNegated())
			{
				conEdges2.add(EdgeDefinition.getVoidEdge());
			}
			else
			{
				conEdges2.add(e);
			}
			if (! e.isSymmetrical())
			{
				// the set of condition-specific edges is only symmetrical when all edges are
				conSymm = false;
			}
		}
		
		if (conEdges.isEmpty())
		{
			conEdges2.add(EdgeDefinition.getVoidEdge());
		}
		
		boolean refNeg = refEdge.isNegated();
		if (refNeg)
			refEdge = EdgeDefinition.getVoidEdge();
		
		diff_edge.makeNegated(false);	// a differential edge is never negated 
		
		// a differential edge is only symmetrical if all input edges are
		boolean refSymm = refEdge.isSymmetrical();
		boolean diffSymm = refSymm && conSymm;
		diff_edge.makeSymmetrical(diffSymm);
		
		//////////// DETERMINE TYPE AND WEIGHT ////////////////
		
		String refCat = getSourceCategory(refEdge.getType());
		
		int countUp = 0;
		int countDown = 0;
		double minDiffWeight = Double.MAX_VALUE;
		double minCumulWeight = Double.MAX_VALUE;
		
		// DEFINE THE COMMON PARENT OF ALL EDGES
		allEdges.addAll(conEdges2);
		allEdges.add(refEdge);
		
		String firstParent = retrieveFirstCommonParent(allEdges, true);
		if (firstParent == null) 	
		{
			return EdgeDefinition.getVoidEdge();
		}
		String firstNeutralParent = firstParent;
		while (firstNeutralParent != null && (posSourceCats.contains(firstNeutralParent) || negSourceCats.contains(firstNeutralParent)))
		{
			firstNeutralParent = retrieveParent(firstNeutralParent);
		}
		
		if (firstNeutralParent == null) 	
		{
			return EdgeDefinition.getVoidEdge();
		}
		
		String baseType = firstNeutralParent;
		boolean unspecified = false;
		
		for (EdgeDefinition conEdge : conEdges2)
		{
			String conCat = getSourceCategory(conEdge.getType());
			double diffWeight = conEdge.getWeight() - refEdge.getWeight();
			double cumWeight = conEdge.getWeight() + refEdge.getWeight();
			
			// refcat is void, concat is not --> increase (unless concat is negative)
			if (refCat.equals(VOID_TYPE) && ! conCat.equals(VOID_TYPE))
			{
				if (negSourceCats.contains(conCat))
				{
					countDown++;
				}
				else
				{
					countUp++;
				}
				diffWeight = conEdge.getWeight();
			}
			
			// refcat is not void, concat is void --> decrease (unless refcat is negative)
			else if (! refCat.equals(VOID_TYPE)  && conCat.equals(VOID_TYPE))
			{
				if (negSourceCats.contains(refCat))
				{
					countUp++;
				}
				else
				{
					countDown++;
				}
				diffWeight = refEdge.getWeight();
			}
			
			// refcat is positive, concat is negative --> decrease
			else if (posSourceCats.contains(refCat) && negSourceCats.contains(conCat))
			{
				countDown++;
				diffWeight = cumWeight;
			}
			
			// refcat is negative, concat is positive --> increase
			else if (negSourceCats.contains(refCat) && posSourceCats.contains(conCat))
			{
				countUp++;
				diffWeight = cumWeight;
			}
			
			else
			{
				if (diffWeight < 0) // decrease
				{
					diffWeight *= -1;
					countDown++;
				}
				else if (diffWeight > 0)	// increase
				{
					countUp++;
				}
				boolean refNeutral = ! (negSourceCats.contains(refCat) || posSourceCats.contains(refCat));
				boolean conNeutral = ! (negSourceCats.contains(conCat) || posSourceCats.contains(conCat));
	
				if ((refNeutral && ! conNeutral) || (conNeutral && ! refNeutral))
				{
					unspecified = true;
				}
			}
			
			minDiffWeight = Math.min(minDiffWeight, diffWeight);
			minCumulWeight = Math.min(minCumulWeight, cumWeight);
		}
		// some are up, some are down -> no general differential edge
		if (countUp > 0 && countDown > 0)
		{
			return EdgeDefinition.getVoidEdge();
		}
		
		// all edges are either all up, or all down
		boolean up = countUp > 0;
		
		String type = "";
		double diffWeight = 0.0;
		
		//type = refCat + "_to_" + conParent;
		//diffWeight = minCumulWeight;
		
		if (up)
		{
			if (diffSymm)
				type = posPrefix_symm;
			else
				type = posPrefix_dir;
		} 
		else
		{
			if (diffSymm)
				type = negPrefix_symm;
			else
				type = negPrefix_dir;
		}
		if (unspecified)
		{
			type += "_unspecified";
		}
		type += "_" + baseType;
		diffWeight = minDiffWeight;
		
		diff_edge.setType(type);
		
		if (diffWeight <= cutoff)
		{
			return EdgeDefinition.getVoidEdge();
		}
		diff_edge.setWeight(diffWeight);
		
		return diff_edge;
	}

	
	////////////// PAINT ALGORITHMS //////////////////////////////////
	
	/**
	 * Assign a specific paint object to a source category (and its children)
	 * 
	 * @param parentCat a category (also representing its children)
	 * @param p the Color object specifying its visual properties
	 * @throws IllegalArgumentException when the either of the arguments are null, when the type was previously assigned to a paint object, 
	 * or when the type is not defined in this ontology
	 */
	protected void addColor(String parentCat, Color p)
	{
		if (parentCat == null || p == null)
		{
			String errormsg = "The provided parent category or the paint object should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (parentSourceCatToColor.containsKey(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') already has a mapped paint object!";
			throw new IllegalArgumentException(errormsg);
		}
		if (! isDefinedSourceCat(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') is not defined in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		parentSourceCatToColor.put(parentCat, p);
	}
	
	/**
	 * Assign a specific arrowhead to a source category (and its children)
	 * 
	 * @param parentCat a category (also representing its children)
	 * @param p the ArrowHead object specifying its visual properties
	 * @throws IllegalArgumentException when the either of the arguments are null, when the type was previously assigned to a paint object, 
	 * or when the type is not defined in this ontology
	 */
	protected void addArrowHead(String parentCat, ArrowHead p)
	{
		if (parentCat == null || p == null)
		{
			String errormsg = "The provided parent category or the ArrowHead object should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (parentSourceCatToArrowHead.containsKey(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') already has a mapped ArrowHead object!";
			throw new IllegalArgumentException(errormsg);
		}
		if (! isDefinedSourceCat(parentCat))
		{
			String errormsg = "The provided parent category ('" + parentCat + "') is not defined in this ontology!";
			throw new IllegalArgumentException(errormsg);
		}
		parentSourceCatToArrowHead.put(parentCat, p);
	}
	
	@Override
	public ArrowHead getDifferentialEdgeArrowHead(String category)
	{
		if (category == null)
		{
			return default_diff_ah;
		}
		if (category.startsWith(negPrefix_symm) || category.startsWith(posPrefix_symm))
		{
			return symm_ah;
		}
		if (category.startsWith(posPrefix_dir))
		{
			return pos_diff_ah;
		}
		if (category.startsWith(negPrefix_dir))
		{
			return neg_diff_ah;
		}
		return default_diff_ah;
	}
	
	
	@Override
	public Color getDifferentialEdgeColor(String category)
	{
		if (category == null)
		{
			return default_diff_paint;
		}
		if (category.startsWith(posPrefix_dir) || category.startsWith(posPrefix_symm))
		{
			return pos_diff_paint;
		}
		if (category.startsWith(negPrefix_symm) || category.startsWith(negPrefix_dir))
		{
			return neg_diff_paint;
		}
		return default_diff_paint;
	}
	
	@Override
	protected ArrowHead getSourceEdgeArrowHead(String edgeType)
	{
		if (isDefinedSourceType(edgeType))
		{
			if (isSymmetricalSourceType(edgeType))
			{
				return symm_ah;
			}
			
			String childCat = getSourceCategory(edgeType);
			ArrowHead foundArrowHead = parentSourceCatToArrowHead.get(edgeType);
			while (foundArrowHead == null && childCat != null)
			{
				String parentCat = retrieveParent(childCat);
				foundArrowHead = parentSourceCatToArrowHead.get(parentCat);
				childCat = parentCat;
			}
			if (foundArrowHead != null)
			{
				return foundArrowHead;
			}
		}
		return neutral_source_ah;
	}
	
	@Override
	protected Color getSourceEdgeColor(String edgeType)
	{
		if (isDefinedSourceType(edgeType))
		{
			String childCat = getSourceCategory(edgeType);
			Color foundColor = parentSourceCatToColor.get(childCat);
			
			while (foundColor == null && childCat != null)
			{
				String parentCat = retrieveParent(childCat);
				foundColor = parentSourceCatToColor.get(parentCat);
				
				childCat = parentCat;
			}
			if (foundColor != null)
			{
				return foundColor;
			}
		}
		return neutral_source_paint;
	}
	
}

package be.svlandeg.diffany.core.algorithms;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.EdgeGenerator;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/**
 * This class does the actual edge-by-edge comparisons, based on the semantic definitions of an ontology, structured as a tree classification system.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeComparison
{

	protected TreeEdgeOntology teo;
	protected EdgeGenerator eg;

	/**
	 * Initialize this object by defining the Edge Ontology which will be used for semantically comparing edges.
	 * 
	 * @param teo the edge ontology defining the interaction semantics
	 */
	public EdgeComparison(TreeEdgeOntology teo)
	{
		this.teo = teo;
		eg = new EdgeGenerator();
	}

	/**
	 * Private class that will be used to store intermediate calculations at all levels of the edge ontology.
	 * Each IntermediateComparison object represents one type in the ontology 
	 * and each possible weight for this type is stored together with the IDs of the supporting networks.
	 */
	protected class IntermediateComparison
	{
		private SortedMap<Double, Set<Integer>> allWeights;
		private String type;

		protected IntermediateComparison(String type)
		{
			allWeights = new TreeMap<Double, Set<Integer>>();
			this.type = type;
		}

		public String toString()
		{
			String result = "intermediate comparison at type " + type + " : support " + getTotalSupport() + " :";
			for (double d : allWeights.keySet())
			{
				result += " " + d + " (" + allWeights.get(d) + ") - ";
			}
			return result;
		}

		public int getTotalSupport()
		{
			int support = 0;
			for (double w : allWeights.keySet())
			{
				support += allWeights.get(w).size();
			}
			return support;
		}
	}

	/**
	 * This method adds one additional piece of support to an intermediate result.
	 * 
	 * @param intermediate the current result
	 * @param networkID the ID of the network which provides this support
	 * @param weight the weight of the edges that further supports this result
	 */
	protected void addResult(IntermediateComparison intermediate, int networkID, double weight)
	{
		Set<Integer> supports = intermediate.allWeights.get(weight);
		if (supports == null)
		{
			supports = new HashSet<Integer>();
		}
		supports.add(networkID);
		intermediate.allWeights.put(weight, supports);
	}

	/**
	 * Add an edge to the set (tree) of intermediate results, by starting with the type of the edge and going up the tree. Each parent gets additional support.
	 * In case the edge is negated, we go down the tree, because the support then travels the other direction.
	 * 
	 * @param e the edge that needs to be added to the intermediate results
	 * @param support the IDs of the networks that provide support for this edge
	 * @param affirmative_results the (flattened) tree of support for non-negated edges
	 * @param negated_results the (flattened) tree of support for negated edges
	 */
	protected void addEdgeToTree(EdgeDefinition e, Set<Integer> support, Map<String, IntermediateComparison> affirmative_results, Map<String, IntermediateComparison> negated_results)
	{
		Set<String> processedCategories = new HashSet<String>(); // make sure we do not process any ontology category twice for the same edge

		if (e.isNegated())
		{
			Set<String> categories = new HashSet<String>();
			categories.add(teo.getSourceCategory(e.getType()));

			while (!categories.isEmpty())
			{
				Set<String> childCategories = new HashSet<String>();

				for (String category : categories)
				{
					if (!processedCategories.contains(category))
					{
						IntermediateComparison intermediate = negated_results.get(category);
						if (intermediate == null)
						{
							intermediate = new IntermediateComparison(category);
						}

						for (int networkID : support)
						{
							addResult(intermediate, networkID, e.getWeight());
						}
						negated_results.put(category, intermediate);

						processedCategories.add(category);
						childCategories.addAll(teo.retrieveCatChildren(category));
					}
				}
				categories = childCategories;
			}
		}
		else
		{
			String category = teo.getSourceCategory(e.getType());
			while (category != null)
			{
				if (!processedCategories.contains(category))
				{
					IntermediateComparison intermediate = affirmative_results.get(category);
					if (intermediate == null)
					{
						intermediate = new IntermediateComparison(category);
					}

					for (int networkID : support)
					{
						addResult(intermediate, networkID, e.getWeight());
					}
					affirmative_results.put(category, intermediate);

					processedCategories.add(category);
					category = teo.retrieveCatParent(category);
				}
				else
				{
					category = null;
				}
			}
		}
	}

	/**
	 * Create all possible edges from an intermediate comparison, by determining the weights which are still supported. Further define the type, symmetry and negation status of the new edge.
	 * If the weight cutoff is not reached, null will be returned
	 */
	private Map<EdgeDefinition, Set<Integer>> createAllEdges(IntermediateComparison inter, boolean final_symm, boolean negation, int overlapNo_cutoff, double weight_cutoff, boolean minOperator)
	{
		Map<EdgeDefinition, Set<Integer>> map = new HashMap<EdgeDefinition, Set<Integer>>();

		// we can take the maximum value, if there is enough total support -> return only one edge with the maximal weight value
		Set<Integer> supports = new HashSet<Integer>();
		if (!minOperator)
		{
			double maxWeight = inter.allWeights.lastKey();
			if (inter.getTotalSupport() >= overlapNo_cutoff && maxWeight >= weight_cutoff)
			{
				for (double w : inter.allWeights.keySet())
				{
					supports.addAll(inter.allWeights.get(w));
				}

				EdgeDefinition overlap_edge = eg.getDefaultEdge();
				overlap_edge.makeSymmetrical(final_symm);
				overlap_edge.makeNegated(negation);
				overlap_edge.setType(inter.type);
				overlap_edge.setWeight(maxWeight);
				map.put(overlap_edge, supports);
				return map;
			}
			return null;
		}

		// When we reach to this point, we need to apply the minimum operator of the edge weights!

		// we need to find all weight values that have enough support (these will be pruned later)
		// starting with the highest weights, their support is passed on to the lower weights. 
		// When the cutoff is reached, corresponding edges will be created

		int accumulatedSupport = 0;

		supports = new HashSet<Integer>();
		TreeSet<Double> allWs = new TreeSet<Double>(inter.allWeights.keySet());
		for (double w : allWs.descendingSet())
		{
			Set<Integer> currentSupport = inter.allWeights.get(w);
			supports.addAll(currentSupport);

			accumulatedSupport += currentSupport.size();
			if (accumulatedSupport >= overlapNo_cutoff && w >= weight_cutoff)
			{
				EdgeDefinition overlap_edge = eg.getDefaultEdge();
				overlap_edge.makeSymmetrical(final_symm);
				overlap_edge.makeNegated(negation);
				overlap_edge.setType(inter.type);
				overlap_edge.setWeight(w);
				map.put(overlap_edge, new HashSet<Integer>(supports));
			}
		}

		return map;
	}


	/**
	 * Method that defines the differential edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns EdgeDefinition.getVoidEdge() when the edge should be deleted (i.e. not present in differential network).
	 * 
	 * An important parameter is overlapNo_cutoff, which determines the amount of support needed for an edge to be included in the consensus condition network. If it equals the number of input networks, all networks need to agree on an edge.
	 * However, if it is smaller, e.g. 3 out of 4, there can be one 'outlier' condition network (potentially a different one for each calculated edge), allowing some noise in the input and creating more robust differential networks.
	 * The overlapNo_cutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting edges (networks).
	 * 
	 * @param refEdge the edge definition in the reference network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param conEdges the edge definitions in the condition-specific networks (can be EdgeDefinition.getVoidEdge() when non-existing)
	 * @param weight_cutoff the minimal weight of a resulting edge for it to be included in the differential network
	 * @param overlapNo_cutoff the minimal number of condition networks (inclusive) that need to have the overlap for it to result to a differential edge
	 * 
	 * @return the edge definition in the differential network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, Collection<EdgeDefinition> conEdges, double weight_cutoff, int overlapNo_cutoff) throws IllegalArgumentException
	{
		EdgeDefinition diff_edge = eg.getDefaultEdge();

		Set<EdgeDefinition> conEdges2 = new HashSet<EdgeDefinition>();
		Set<EdgeDefinition> allEdges = new HashSet<EdgeDefinition>();

		Set<String> posSourceCats = teo.getAllPosSourceCategories();
		Set<String> negSourceCats = teo.getAllNegSourceCategories();

		//////////// DEAL WITH SYMMETRY AND NEGATION ////////////////
		boolean conSymm = true;
		for (EdgeDefinition e : conEdges)
		{
			if (!e.isSymmetrical())
			{
				// the set of condition-specific edges is only symmetrical when all edges are
				conSymm = false;
			}
		}
		for (EdgeDefinition e : conEdges)
		{
			// negated edges are set to void
			if (e.isNegated())
			{
				conEdges2.add(eg.getVoidEdge(conSymm));
			}
			else
			{
				conEdges2.add(e);
			}
		}

		if (conEdges.isEmpty())
		{
			conEdges2.add(eg.getVoidEdge(conSymm));
		}

		boolean refNeg = refEdge.isNegated();
		if (refNeg)
		{
			refEdge = eg.getVoidEdge(conSymm);
		}

		diff_edge.makeNegated(false); // a differential edge is never negated 

		// a differential edge is only symmetrical if all input edges are
		boolean refSymm = refEdge.isSymmetrical();
		boolean diffSymm = refSymm && conSymm;
		diff_edge.makeSymmetrical(diffSymm);

		//////////// DETERMINE TYPE AND WEIGHT ////////////////

		String refCat = teo.getSourceCategory(refEdge.getType());

		int countUp = 0;
		int countDown = 0;
		double minDiffWeight = Double.MAX_VALUE;
		double minCumulWeight = Double.MAX_VALUE;

		// DEFINE THE COMMON PARENT OF ALL EDGES
		allEdges.addAll(conEdges2);
		allEdges.add(refEdge);

		Set<String> cats = new HashSet<String>();
		int countEmpty = 0;

		for (EdgeDefinition e : allEdges)
		{
			String cat = teo.getSourceCategory(e.getType());
			if (cat.equals(teo.getVoidCategory(e.isSymmetrical())))
			{
				countEmpty++;
			}
			else
			{
				cats.add(cat);
			}
		}

		String firstParent = null;
		if (countEmpty != allEdges.size())
		{
			firstParent = teo.retrieveFirstCommonParent(cats);
		}

		if (firstParent == null)
		{
			return eg.getVoidEdge(conSymm);
		}
		String firstNeutralParent = firstParent;
		while (firstNeutralParent != null && (posSourceCats.contains(firstNeutralParent) || negSourceCats.contains(firstNeutralParent)))
		{
			firstNeutralParent = teo.retrieveCatParent(firstNeutralParent);
		}

		if (firstNeutralParent == null)
		{
			return eg.getVoidEdge(conSymm);
		}

		String baseType = firstNeutralParent;
		boolean unspecified = false;

		String VOID_CAT = teo.getVoidCategory(conSymm);

		for (EdgeDefinition conEdge : conEdges2)
		{
			String conCat = teo.getSourceCategory(conEdge.getType());
			double diffWeight = conEdge.getWeight() - refEdge.getWeight();
			double cumWeight = conEdge.getWeight() + refEdge.getWeight();

			// refcat is void, concat is not --> increase (unless concat is negative)
			if (refCat.equals(VOID_CAT) && !conCat.equals(VOID_CAT))
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
			else if (!refCat.equals(VOID_CAT) && conCat.equals(VOID_CAT))
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
				else if (diffWeight > 0) // increase
				{
					countUp++;
				}
				boolean refNeutral = !(negSourceCats.contains(refCat) || posSourceCats.contains(refCat));
				boolean conNeutral = !(negSourceCats.contains(conCat) || posSourceCats.contains(conCat));

				if ((refNeutral && !conNeutral) || (conNeutral && !refNeutral))
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
			return eg.getVoidEdge(conSymm);
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
				type = teo.getPosPrefix_symm();
			else
				type = teo.getPosPrefix_dir();
		}
		else
		{
			if (diffSymm)
				type = teo.getNegPrefix_symm();
			else
				type = teo.getNegPrefix_dir();
		}
		if (unspecified)
		{
			type += teo.getUnspecifiedPrefix();
		}
		type += baseType;
		diffWeight = minDiffWeight;

		diff_edge.setType(type);

		if (diffWeight <= weight_cutoff)
		{
			return eg.getVoidEdge(conSymm);
		}
		diff_edge.setWeight(diffWeight);
		return diff_edge;
	}

	/**
	 * Method that defines the overlapping edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns an empty set when the edge should be deleted (i.e. not present in the overlapping network).
	 * 
	 * An important parameter is overlapNo_cutoff, which determines the amount of support needed for an edge to be included in the overlap network. If it equals the number of input networks, all networks need to agree on an edge.
	 * However, if it is smaller, e.g. 3 out of 4, there can be one 'outlier' network (potentially a different one for each calculated edge), allowing some noise in the input and creating more robust overlap networks.
	 * The overlapNo_cutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting edges (networks).
	 * 
	 * @param edges the original edge definitions (should not be null or empty!)
	 * @param supports the network IDs which support the corresponding edges in the other list 
	 * @param overlapNo_cutoff the minimal number of networks (inclusive) that need to have the overlap for it to be included
	 * @param weight_cutoff the minimal value of a resulting edge for it to be included in the overlapping network
	 * @param minOperator whether or not to take the minimum of the edge weights - if false, the maximum is taken
	 * 
	 * @return the edge definitions in the overlapping network, or an empty set, but never null
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public Map<EdgeDefinition, Set<Integer>> getOverlapEdge(List<EdgeDefinition> edges, List<Set<Integer>> supports, int overlapNo_cutoff, double weight_cutoff, boolean minOperator) throws IllegalArgumentException
	{
		Map<EdgeDefinition, Set<Integer>> overlaps = new HashMap<EdgeDefinition, Set<Integer>>();

		if (edges == null || edges.isEmpty())
		{
			String errormsg = "The list of edges should not be null or empty!";
			throw new IllegalArgumentException(errormsg);
		}
		if (supports == null || supports.isEmpty() || supports.size() != edges.size())
		{
			String errormsg = "The list of supporting networks should be as large as the edge list!";
			throw new IllegalArgumentException(errormsg);
		}
		int countEdges = edges.size();

		// 0. CHECK SYMMETRY //
		int countSymmetrical = 0;

		for (EdgeDefinition e : edges)
		{
			if (e.isSymmetrical())
			{
				countSymmetrical++;
			}
		}

		if (countSymmetrical != countEdges && countSymmetrical != 0)
		{
			String errormsg = "The set of edges should either be all symmetrical, or all asymmetrical - clean the input first with NetworkCleaning.unifyDirection !";
			throw new IllegalArgumentException(errormsg);
		}

		boolean final_symm = countSymmetrical == countEdges;

		// 1. PROCESS ALL EDGES ONE BY ONE AND ADD THEM TO THE INTERMEDIATE RESULTS //

		Map<String, IntermediateComparison> affirmative_results = new HashMap<String, IntermediateComparison>();
		Map<String, IntermediateComparison> negated_results = new HashMap<String, IntermediateComparison>();

		for (int i = 0; i < edges.size(); i++)
		{
			EdgeDefinition e = edges.get(i);
			Set<Integer> support = supports.get(i);
			addEdgeToTree(e, support, affirmative_results, negated_results);
		}

		// 2. GO THROUGH THE WHOLE ONTOLOGY TREE AND COLLECT THE RESULTS, USING THE OVERLAP NO CUTOFF PARAMETER

		for (String cat : teo.getAllSourceCategories())
		{
			IntermediateComparison aff_result = affirmative_results.get(cat);
			boolean aff = false;
			if (aff_result != null && aff_result.getTotalSupport() >= overlapNo_cutoff)
			{
				aff = true;
			}

			IntermediateComparison neg_result = negated_results.get(cat);
			boolean neg = false;
			if (neg_result != null && neg_result.getTotalSupport() >= overlapNo_cutoff)
			{
				neg = true;
			}

			if (aff)
			{
				Map<EdgeDefinition, Set<Integer>> map = createAllEdges(aff_result, final_symm, false, overlapNo_cutoff, weight_cutoff, minOperator);
				if (map != null)
				{
					for (EdgeDefinition overlap_edge : map.keySet())
					{
						Set<Integer> theseSupports = map.get(overlap_edge);
						overlaps.put(overlap_edge, theseSupports);
					}
				}
			}
			if (neg)
			{
				Map<EdgeDefinition, Set<Integer>> map = createAllEdges(neg_result, final_symm, true, overlapNo_cutoff, weight_cutoff, minOperator);
				if (map != null)
				{
					for (EdgeDefinition overlap_edge : map.keySet())
					{
						Set<Integer> theseSupports = map.get(overlap_edge);
						overlaps.put(overlap_edge, theseSupports);
					}
				}
			}
		}

		return overlaps;
	}

}

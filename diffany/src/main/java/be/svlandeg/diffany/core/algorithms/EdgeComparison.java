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
	 */
	protected class IntermediateComparison
	{
		private SortedMap<Double, Integer> allWeights;
		private String type;

		private int support; // the number of edges matching this level of the ontology

		protected IntermediateComparison(String type)
		{
			allWeights = new TreeMap<Double, Integer>();
			support = 0;
			this.type = type;
		}

		public String toString()
		{
			String result = "intermediate comparison at type " + type + " : support " + support + " :";
			for (double d : allWeights.keySet())
			{
				result += " " + d + " (" + allWeights.get(d) + ") - ";
			}
			return result;
		}
	}

	/**
	 * This method adds one additional piece of support to an intermediate result.
	 * 
	 * @param intermediate the current result
	 * @param weight the weight of the edges that further supports this result
	 */
	protected void addResult(IntermediateComparison intermediate, double weight)
	{
		Integer count = intermediate.allWeights.get(weight);
		if (count == null)
		{
			count = 0;
		}
		count++;
		intermediate.allWeights.put(weight, count);
		intermediate.support = intermediate.support + 1;
	}

	/**
	 * Determine the optimal weight to be given to an edge, depending on the needed level of support and the min/max operator
	 * Return null when an appropriate weight could not be found
	 * 
	 * @param weight_cutoff the minimal value of a resulting edge for it to be included in the overlapping network
	 * 
	 * TODO javadoc
	 */
	protected Double determineWeight(IntermediateComparison intermediate, int overlapNo_cutoff, double weight_cutoff, boolean minOperator)
	{
		// we can take the maximum value, if there is enough support
		if (!minOperator)
		{
			double maxWeight = intermediate.allWeights.lastKey();
			if (intermediate.support >= overlapNo_cutoff && maxWeight >= weight_cutoff)
			{
				return maxWeight;
			}
			return null;
		}

		// we need to find the maximal minimum value that has enough support
		// starting with the highest weights, their support is passed on to the lower weights until the cutoff is reached
		int accumulatedSupport = 0;

		TreeSet<Double> allWs = new TreeSet<Double>(intermediate.allWeights.keySet());
		for (double w : allWs.descendingSet())
		{
			int currentSupport = intermediate.allWeights.get(w);
			accumulatedSupport += currentSupport;
			if (accumulatedSupport >= overlapNo_cutoff && w >= weight_cutoff)
			{
				return w; // the first (highest) value that has enough support
			}
		}

		return null;
	}

	/**
	 * Add an edge to the set (tree) of intermediate results, by starting with the type of the edge and going up the tree. Each parent gets additional support.
	 * In case the edge is negated, we go down the tree, because the support then travels the other direction.
	 * 
	 * @param e the edge that needs to be added to the intermediate results
	 * @param affirmative_results the (flattened) tree of support for non-negated edges
	 * @param negated_results the (flattened) tree of support for negated edges
	 */
	protected void addEdgeToTree(EdgeDefinition e, Map<String, IntermediateComparison> affirmative_results, Map<String, IntermediateComparison> negated_results)
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

						addResult(intermediate, e.getWeight());
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

					addResult(intermediate, e.getWeight());
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
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, Collection<EdgeDefinition> conEdges, double cutoff) throws IllegalArgumentException
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

		if (diffWeight <= cutoff)
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
		if (supports == null || supports.isEmpty())
		{
			String errormsg = "The list of supporting networks should not be null or empty!";
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

		for (EdgeDefinition e : edges)
		{
			addEdgeToTree(e, affirmative_results, negated_results);
		}

		// 2. GO THROUGH THE WHOLE ONTOLOGY TREE AND COLLECT THE RESULTS, USING THE OVERLAP NO CUTOFF PARAMETER

		for (String cat : teo.getAllSourceCategories())
		{
			IntermediateComparison aff_result = affirmative_results.get(cat);
			boolean aff = false;
			if (aff_result != null && aff_result.support >= overlapNo_cutoff)
			{
				aff = true;
			}
			
			IntermediateComparison neg_result = negated_results.get(cat);
			boolean neg = false;
			if (neg_result != null && neg_result.support >= overlapNo_cutoff)
			{
				neg = true;
			}
			
			if (aff)
			{
				EdgeDefinition overlap_edge = createEdge(aff_result, final_symm, false, overlapNo_cutoff, weight_cutoff, minOperator);
				if (overlap_edge != null)
				{
					overlaps.put(overlap_edge, null);	// TODO
				}
			}
			if (neg)
			{
				EdgeDefinition overlap_edge = createEdge(neg_result, final_symm, true, overlapNo_cutoff, weight_cutoff, minOperator);
				if (overlap_edge != null)
				{
					overlaps.put(overlap_edge, null);	// TODO
				}
			}
		}

		return overlaps;
	}
	
	/**
	 * TODO
	 * 
	 */
	private EdgeDefinition createEdge(IntermediateComparison inter, boolean final_symm, boolean negation, int overlapNo_cutoff, double weight_cutoff, boolean minOperator)
	{
		EdgeDefinition overlap_edge = eg.getDefaultEdge();
		overlap_edge.makeSymmetrical(final_symm);
		overlap_edge.makeNegated(negation);
		overlap_edge.setType(inter.type);
		Double weight = determineWeight(inter, overlapNo_cutoff, weight_cutoff, minOperator);
		if (weight != null)
		{
			overlap_edge.setWeight(weight);
			return overlap_edge;
		}
		return null;
	}

}

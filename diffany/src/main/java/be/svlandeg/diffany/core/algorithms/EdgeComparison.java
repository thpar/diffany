package be.svlandeg.diffany.core.algorithms;

import java.util.ArrayList;
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
	 * @param negated_results the (flattened) tree of support for negated edges (can be null if the edge is not negated)
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
	 * Create all possible edges from an intermediate comparison, by determining the weights which are still supported. 
	 * Further define the type, symmetry and negation status of the new edge.
	 * If the weight cutoff is not reached, null will be returned
	 */
	private Map<EdgeDefinition, Set<Integer>> createAllEdges(IntermediateComparison inter, boolean final_symm, boolean negation, int overlapNo_cutoff, double weight_min, double weight_max, boolean minOperator)
	{
		Map<EdgeDefinition, Set<Integer>> map = new HashMap<EdgeDefinition, Set<Integer>>();

		// we can take the maximum value, if there is enough total support -> return only one edge with the maximal weight value
		Set<Integer> supports = new HashSet<Integer>();
		if (!minOperator)
		{
			double maxWeight = inter.allWeights.lastKey();
			if (inter.getTotalSupport() >= overlapNo_cutoff && maxWeight >= weight_min && maxWeight < weight_max)
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
			if (accumulatedSupport >= overlapNo_cutoff && w >= weight_min && w < weight_max)
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
	 * @param supports the network IDs which support the corresponding conEdges
	 * @param weight_cutoff the minimal weight of a resulting edge for it to be included in the differential network
	 * @param overlapNo_cutoff the minimal number of condition networks (inclusive) that need to have the overlap for it to result to a differential edge
	 * 
	 * @return the edge definition in the differential network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, List<EdgeDefinition> conEdges, List<Set<Integer>> supports, int overlapNo_cutoff, double weight_cutoff) throws IllegalArgumentException
	{
		if (conEdges == null || conEdges.isEmpty())
		{
			String errormsg = "The list of conditional edges should not be null or empty!";
			throw new IllegalArgumentException(errormsg);
		}
		if (supports == null || supports.isEmpty() || supports.size() != conEdges.size())
		{
			String errormsg = "The list of supporting conditional networks should be as large as the conditional edge list!";
			throw new IllegalArgumentException(errormsg);
		}

		EdgeDefinition diff_edge = eg.getDefaultEdge();

		int originalEdges = conEdges.size() + 1;

		// 0. DEAL WITH SYMMETRY AND NEGATION //

		// create the set of condition edges by removing the negated ones, and count symmetry while we're at it
		List<EdgeDefinition> conEdgesPos = new ArrayList<EdgeDefinition>();
		int countSymmetrical = 0;

		for (EdgeDefinition e : conEdges)
		{
			if (e.isSymmetrical())
			{
				countSymmetrical++;
			}
			// negated edges are set to void
			if (!e.isNegated())
			{
				conEdgesPos.add(e);
			}
			else
			{
				conEdgesPos.add(eg.getVoidEdge(e.isSymmetrical()));
			}
		}
		if (refEdge.isSymmetrical())
		{
			countSymmetrical++;
		}

		boolean refNeg = refEdge.isNegated();
		if (refNeg)
		{
			refEdge = eg.getVoidEdge(refEdge.isSymmetrical());
		}

		if (countSymmetrical != originalEdges && countSymmetrical != 0)
		{
			String errormsg = "The set of reference and condition edges should either be all symmetrical, or all asymmetrical - clean the input first with NetworkCleaning.unifyDirection !";
			throw new IllegalArgumentException(errormsg);
		}
		// a differential edge is only symmetrical if all input edges are
		boolean final_symm = countSymmetrical == originalEdges;
		diff_edge.makeSymmetrical(final_symm);

		// a differential edge is never negated 
		diff_edge.makeNegated(false);

		// 1. PROCESS ALL EDGES ONE BY ONE AND ADD THEM TO THE INTERMEDIATE RESULTS //

		Map<String, IntermediateComparison> con_results = new HashMap<String, IntermediateComparison>();

		for (int i = 0; i < conEdgesPos.size(); i++)
		{
			EdgeDefinition e = conEdgesPos.get(i);
			Set<Integer> support = supports.get(i);
			addEdgeToTree(e, support, con_results, null);
		}

		// 2. GO THROUGH THE WHOLE ONTOLOGY TREE AND COLLECT THE CONSENSUS OVERLAP RESULTS  //

		Set<EdgeDefinition> overlaps = new HashSet<EdgeDefinition>();
		double refWeight = refEdge.getWeight();
		boolean minOperator = true;

		for (String cat : teo.getAllSourceCategories())
		{
			IntermediateComparison con_result = con_results.get(cat);

			if (con_result != null && con_result.getTotalSupport() >= overlapNo_cutoff)
			{
				Map<EdgeDefinition, Set<Integer>> map_below = createAllEdges(con_result, final_symm, false, overlapNo_cutoff, Double.NEGATIVE_INFINITY, refWeight, minOperator);
				Map<EdgeDefinition, Set<Integer>> map_above = createAllEdges(con_result, final_symm, false, overlapNo_cutoff, refWeight, Double.POSITIVE_INFINITY, minOperator);

				if (!map_below.isEmpty() && !map_above.isEmpty())
				{
					// there can be no differential edge because there is evidence for both higher and lower weights
					return eg.getVoidEdge(final_symm);
				}
				// keep only the consensus condition edge with the highest weight
				else if (!map_below.isEmpty())
				{
					double max = Double.NEGATIVE_INFINITY;
					for (EdgeDefinition overlap_edge : map_below.keySet())
					{
						max = Math.max(max, overlap_edge.getWeight());
					}
					for (EdgeDefinition overlap_edge : map_below.keySet())
					{
						if (overlap_edge.getWeight() == max)
						{
							overlaps.add(overlap_edge);
						}
					}
				}
				// keep only the consensus condition edge with the highest weight
				else if (!map_above.isEmpty())
				{
					double max = Double.NEGATIVE_INFINITY;
					for (EdgeDefinition overlap_edge : map_above.keySet())
					{
						max = Math.max(max, overlap_edge.getWeight());
					}
					for (EdgeDefinition overlap_edge : map_above.keySet())
					{
						if (overlap_edge.getWeight() == max)
						{
							overlaps.add(overlap_edge);
						}
					}
				}
			}
		}

		if (overlaps.size() == 0)
		{
			overlaps.add(eg.getVoidEdge(final_symm));
		}
		while (overlaps.size() > 1)
		{
			Set<EdgeDefinition> newOverlaps = new HashSet<EdgeDefinition>();
			newOverlaps.addAll(overlaps);
			EdgeDefinition consensusConEdgeChild = overlaps.iterator().next();
			EdgeDefinition consensusConEdgeParent = overlaps.iterator().next();
			if (teo.isSourceCatChildOf(consensusConEdgeChild.getType(), consensusConEdgeParent.getType()) > -1)
			{
				newOverlaps.remove(consensusConEdgeParent);
			}
			overlaps = newOverlaps;
		}

		EdgeDefinition consensusConEdge = overlaps.iterator().next();
		System.out.println("consensusConEdge " + consensusConEdge);

		// 3. DETERMINE THE FINAL TYPE AND WEIGHT //

		String refCat = teo.getSourceCategory(refEdge.getType());

		boolean isDown = false;
		

		// DEFINE THE COMMON PARENT OF ALL EDGES
		Set<EdgeDefinition> allEdges = new HashSet<EdgeDefinition>();
		allEdges.add(consensusConEdge);
		allEdges.add(refEdge);

		Set<String> posSourceCats = teo.getAllPosSourceCategories();
		Set<String> negSourceCats = teo.getAllNegSourceCategories();

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
			return eg.getVoidEdge(final_symm);
		}
		String firstNeutralParent = firstParent;
		while (firstNeutralParent != null && (posSourceCats.contains(firstNeutralParent) || negSourceCats.contains(firstNeutralParent)))
		{
			firstNeutralParent = teo.retrieveCatParent(firstNeutralParent);
		}

		if (firstNeutralParent == null)
		{
			return eg.getVoidEdge(final_symm);
		}

		String baseType = firstNeutralParent;
		boolean unspecified = false;

		String VOID_CAT = teo.getVoidCategory(final_symm);

		String conCat = teo.getSourceCategory(consensusConEdge.getType());
		
		double weightDifference = consensusConEdge.getWeight() - refEdge.getWeight();
		double weightSum = consensusConEdge.getWeight() + refEdge.getWeight();
		double finalDiffWeight = 0;

		// refcat is void, concat is not --> increase (unless concat is negative)
		if (refCat.equals(VOID_CAT) && !conCat.equals(VOID_CAT))
		{
			isDown = negSourceCats.contains(conCat);
			finalDiffWeight = consensusConEdge.getWeight();
		}

		// refcat is not void, concat is void --> decrease (unless refcat is negative)
		else if (!refCat.equals(VOID_CAT) && conCat.equals(VOID_CAT))
		{
			isDown = !negSourceCats.contains(refCat);
			finalDiffWeight = refEdge.getWeight();
		}

		// refcat is positive, concat is negative --> decrease
		else if (posSourceCats.contains(refCat) && negSourceCats.contains(conCat))
		{
			isDown = true;
			finalDiffWeight = weightSum;
		}

		// refcat is negative, concat is positive --> increase
		else if (negSourceCats.contains(refCat) && posSourceCats.contains(conCat))
		{
			isDown = false;
			finalDiffWeight = weightSum;
		}

		else
		{
			finalDiffWeight = weightDifference;
			if (finalDiffWeight < 0) // decrease
			{
				finalDiffWeight *= -1;
				isDown = true;
			}
			else if (finalDiffWeight > 0) // increase
			{
				isDown = false;
			}
			boolean refNeutral = !(negSourceCats.contains(refCat) || posSourceCats.contains(refCat));
			boolean conNeutral = !(negSourceCats.contains(conCat) || posSourceCats.contains(conCat));

			if ((refNeutral && !conNeutral) || (conNeutral && !refNeutral))
			{
				unspecified = true;
			}
		}
		

		String type = "";

		if (isDown)
		{
			if (final_symm)
			{
				type = teo.getNegPrefix_symm();
			}
			else
			{
				type = teo.getNegPrefix_dir();
			}
		}
		else
		{
			if (final_symm)
			{
				type = teo.getPosPrefix_symm();
			}
			else
			{
				type = teo.getPosPrefix_dir();
			}
		}
		if (unspecified)
		{
			type += teo.getUnspecifiedPrefix();
		}
		type += baseType;

		diff_edge.setType(type);

		if (finalDiffWeight <= weight_cutoff)
		{
			return eg.getVoidEdge(final_symm);
		}
		diff_edge.setWeight(finalDiffWeight);
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

		// 2. GO THROUGH THE WHOLE ONTOLOGY TREE AND COLLECT THE RESULTS, USING THE OVERLAP NO CUTOFF PARAMETER //

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
				Map<EdgeDefinition, Set<Integer>> map = createAllEdges(aff_result, final_symm, false, overlapNo_cutoff, weight_cutoff, Double.POSITIVE_INFINITY, minOperator);
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
				Map<EdgeDefinition, Set<Integer>> map = createAllEdges(neg_result, final_symm, true, overlapNo_cutoff, weight_cutoff, Double.POSITIVE_INFINITY, minOperator);
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

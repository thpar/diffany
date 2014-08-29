package be.svlandeg.diffany.core.algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
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

		// generic, zero-weight (non-existing) edges give support 0 to all categories. Zero weight edges are supposedly 'affirmative'.
		if (e.getWeight() == 0)
		{
			Boolean symm_generic = null;
			if (e.getType().equals(new EdgeGenerator().getVoidEdge(true).getType()))
			{
				symm_generic = true;
			}
			if (e.getType().equals(new EdgeGenerator().getVoidEdge(false).getType()))
			{
				symm_generic = false;
			}
			if (symm_generic != null)
			{
				// TODO this should be limited to the relevant subtree!
				Set<String> categories = teo.getAllSourceCategories(true);
				for (String category : categories)
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
					}
				}
			}

		}

		// negated edges travel down the tree, e.g. 'no ptm' also means 'no phosphorylation'
		else if (e.isNegated())
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
		// affirmative edges travel up the tree, e.g. phosphorylation also mean ptm
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
	 * The minimum and maximum weights are *exclusive*.
	 * 
	 * Further define the type, symmetry and negation status of the new edge.
	 * If the weight cutoff is not reached, an empty map will be returned
	 */
	private Map<EdgeDefinition, Set<Integer>> createAllEdges(IntermediateComparison inter, boolean final_symm, boolean negation, int supportingCutoff, double weight_min, double weight_max, boolean minOperator)
	{
		Map<EdgeDefinition, Set<Integer>> map = new HashMap<EdgeDefinition, Set<Integer>>();
		
		int accumulatedSupport = 0;
		Set<Integer> supports = new HashSet<Integer>();
		TreeSet<Double> allWs = new TreeSet<Double>(inter.allWeights.keySet());
		
		// We start from the highest possible weights and keep going down until enough accumulated support is found
		// For the maximum operator, the highest weight is chosen that fits into the provided weight interval
		
		if (!minOperator)
		{
			double highestW = Double.NEGATIVE_INFINITY;
			for (double w : allWs.descendingSet())
			{
				if (w > weight_min && w < weight_max)
				{
					highestW = Math.max(highestW, w);
					Set<Integer> currentSupport = inter.allWeights.get(w);
					supports.addAll(currentSupport);
					accumulatedSupport += currentSupport.size();

					if ((accumulatedSupport >= supportingCutoff))
					{
						EdgeDefinition consensusEdge = eg.getDefaultEdge();
						consensusEdge.makeSymmetrical(final_symm);
						consensusEdge.makeNegated(negation);
						consensusEdge.setType(inter.type);
						consensusEdge.setWeight(highestW);
						map.put(consensusEdge, new HashSet<Integer>(supports));
					}
				}
			}
			return map;
		}

		// For the minimum operator, the first lower weight is chosen that has enough accumulated support
		
		for (double w : allWs.descendingSet())
		{
			if (w > weight_min && w < weight_max)
			{
				Set<Integer> currentSupport = inter.allWeights.get(w);
				supports.addAll(currentSupport);
				accumulatedSupport += currentSupport.size();

				if ((accumulatedSupport >= supportingCutoff))
				{
					EdgeDefinition consensusEdge = eg.getDefaultEdge();
					consensusEdge.makeSymmetrical(final_symm);
					consensusEdge.makeNegated(negation);
					consensusEdge.setType(inter.type);
					consensusEdge.setWeight(w);
					map.put(consensusEdge, new HashSet<Integer>(supports));
				}
			}
		}

		return map;
	}
	
	/**
	 * Create all possible edges from an intermediate comparison, by determining the exact weight we want the supporting edge for.
	 * 
	 * Further define the type, symmetry and negation status of the new edge
	 */
	private Map<EdgeDefinition, Set<Integer>> createAllEdges(IntermediateComparison inter, boolean final_symm, boolean negation, int supportingCutoff, double weight)
	{
		Map<EdgeDefinition, Set<Integer>> map = new HashMap<EdgeDefinition, Set<Integer>>();
				
		Set<Integer> currentSupport = inter.allWeights.get(weight);
		int size = 0;
		if (currentSupport != null)
		{
			size = currentSupport.size();
		}

		if (size >= supportingCutoff)
		{
			EdgeDefinition consensusEdge = eg.getDefaultEdge();
			consensusEdge.makeSymmetrical(final_symm);
			consensusEdge.makeNegated(negation);
			consensusEdge.setType(inter.type);
			consensusEdge.setWeight(weight);
			map.put(consensusEdge, new HashSet<Integer>(currentSupport));
		}
		return map;
	}

	/**
	 * Method that defines the differential edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns EdgeDefinition.getVoidEdge() when the edge should be deleted (i.e. not present in differential network).
	 * 
	 * An important parameter is supportingCutoff, which determines the amount of support needed for an edge to be included in the consensus condition network. If it equals the number of input networks, all networks need to agree on an edge.
	 * However, if it is smaller, e.g. 3 out of 4, there can be one 'outlier' condition network (potentially a different one for each calculated edge), allowing some noise in the input and creating more robust differential networks.
	 * The supportingCutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting edges (networks).
	 * 
	 * @param refEdge the edge definition in the reference network (can be a EdgeDefinition.getVoidEdge() when non-existing)
	 * @param conEdges the edge definitions in the condition-specific networks (can be EdgeDefinition.getVoidEdge() when non-existing)
	 * @param supports the network IDs which support the corresponding conEdges
	 * @param weightCutoff the minimal weight of a resulting edge for it to be included in the differential network
	 * @param supportingCutoff the minimal number of condition networks (inclusive) that need to have the consensus for it to result to a differential edge
	 * 
	 * @return the edge definition in the differential network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, List<EdgeDefinition> conEdges, List<Set<Integer>> supports, int supportingCutoff, double weightCutoff) throws IllegalArgumentException
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

		if (refEdge.isNegated() || refEdge.getWeight() == 0)
		{
			refEdge = eg.getVoidEdge(refEdge.isSymmetrical());
		}
		double refWeight = refEdge.getWeight();

		boolean allEmpty = true;
		if (refWeight != 0)
		{
			allEmpty = false;
		}

		EdgeDefinition diff_edge = eg.getDefaultEdge();

		int originalEdges = conEdges.size() + 1;

		// 0. DEAL WITH SYMMETRY AND NEGATION //

		// create the set of condition edges by removing the negated ones, and count symmetry while we're at it
		List<EdgeDefinition> conEdgesPos = new ArrayList<EdgeDefinition>();
		int countSymmetrical = 0;

		for (EdgeDefinition e : conEdges)
		{
			if (e.getWeight() != 0)
			{
				allEmpty = false;
			}
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

		if (countSymmetrical != originalEdges && countSymmetrical != 0)
		{
			String errormsg = "The set of reference and condition edges should either be all symmetrical, or all asymmetrical - clean the input first with NetworkCleaning.unifyDirection !";
			throw new IllegalArgumentException(errormsg);
		}
		// a differential edge is only symmetrical if all input edges are
		boolean final_symm = countSymmetrical == originalEdges;
		diff_edge.makeSymmetrical(final_symm);
		if (allEmpty)
		{
			return eg.getVoidEdge(final_symm);
		}

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

		// 2. GO THROUGH THE WHOLE ONTOLOGY TREE AND COLLECT THE CONSENSUS RESULTS  //

		Set<EdgeDefinition> consensusAll = new HashSet<EdgeDefinition>();
		Set<EdgeDefinition> consensusClean = new HashSet<EdgeDefinition>();

		for (String cat : teo.getAllSourceCategories(true))
		{
			IntermediateComparison con_result = con_results.get(cat);
			if (con_result != null && con_result.getTotalSupport() >= supportingCutoff)
			{
				// all edges with weight below the reference weight: take the maximum of the condition edges to determine the minimal consensus decrease
				Map<EdgeDefinition, Set<Integer>> map_below = createAllEdges(con_result, final_symm, false, supportingCutoff, Double.NEGATIVE_INFINITY, refWeight, false);
				
				// all edges with weight above the reference weight: take the minimum of the condition edges to determine the minimal consensus increase
				Map<EdgeDefinition, Set<Integer>> map_above = createAllEdges(con_result, final_symm, false, supportingCutoff, refWeight, Double.POSITIVE_INFINITY, true);
				
				Map<EdgeDefinition, Set<Integer>> map_same = createAllEdges(con_result, final_symm, false, supportingCutoff, refWeight);
				
				if (map_below.isEmpty() && map_above.isEmpty())
				{
					if (map_same.isEmpty())
					{
					// there can be no differential edge because the evidence is spread between higher and lower weights
						return eg.getVoidEdge(final_symm);
					}
					if (! map_same.isEmpty() && refWeight != 0)
					{
						// there can be no differential edge because the conditional weights are equal to the ref Edge
						// in case the edge weight was 0, we first continue the search in other categories
						return eg.getVoidEdge(final_symm);
					}
				}
				else if (!map_below.isEmpty() && !map_above.isEmpty())
				{
					// there can be no differential edge because there is evidence for both higher and lower weights
					return eg.getVoidEdge(final_symm);
				}
				// keep only the consensus condition edge with the highest weight (below threshold)
				else if (!map_below.isEmpty())
				{
					for (EdgeDefinition consensusEedge : map_below.keySet())
					{
						consensusAll.add(consensusEedge);
					}
				}
				// keep only the consensus condition edge with the highest weight (above threshold)
				else if (!map_above.isEmpty())
				{
					for (EdgeDefinition consensusEedge : map_above.keySet())
					{
						consensusAll.add(consensusEedge);
					}
				}
			}
		}

		double max = Double.NEGATIVE_INFINITY;
		for (EdgeDefinition consensusEedge : consensusAll)
		{
			max = Math.max(max, consensusEedge.getWeight());
		}
		if (max > 0)
		{
			for (EdgeDefinition consensusEedge : consensusAll)
			{
				if (consensusEedge.getWeight() == max)
				{
					consensusClean.add(consensusEedge);
				}
			}
		}

		if (consensusClean.size() == 0)
		{
			consensusClean.add(eg.getVoidEdge(final_symm));
		}
		// take the most specific condition category
		while (consensusClean.size() > 1)
		{
			// TODO: this sometimes generates an infinite loop!!!

			Iterator<EdgeDefinition> it = consensusClean.iterator();
			EdgeDefinition consensusEdge0 = it.next();
			EdgeDefinition consensusEdge1 = it.next();

			// we copy all edges, except for the possible "child" or "parent"
			Set<EdgeDefinition> newEdges = new HashSet<EdgeDefinition>();
			while (it.hasNext())
			{
				newEdges.add(it.next());
			}

			if (teo.isSourceCatChildOf(consensusEdge0.getType(), consensusEdge1.getType()) >= 0)
			{
				// we only add the child (or equal) to the new list
				newEdges.add(consensusEdge0);
			}
			else if (teo.isSourceCatChildOf(consensusEdge1.getType(), consensusEdge0.getType()) > 0)
			{
				// we only add the child to the new list
				newEdges.add(consensusEdge1);
			}
			else
			{
				newEdges.add(consensusEdge0);
				newEdges.add(consensusEdge1);
			}
			consensusClean = newEdges;
		}

		EdgeDefinition consensusConEdge = consensusClean.iterator().next();

		double conWeight = consensusConEdge.getWeight();

		// 3. DETERMINE THE FINAL TYPE AND WEIGHT //

		String refCat = teo.getSourceCategory(refEdge.getType());

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
			if (e.getWeight() == 0)
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

		Boolean direction = null; // null=unspecified, true=increase, false=decrease
		boolean unspecified = false;

		String conCat = teo.getSourceCategory(consensusConEdge.getType());

		double finalDiffWeight = 0;

		// We assume that void edges are of weight 0 !

		// both edges are void: return an empty differential edge
		if (refWeight == 0 && conWeight == 0)
		{
			return eg.getVoidEdge(final_symm);
		}

		// ref edge is void, con edge is not --> differential increase (unless concat is negative)
		else if (refWeight == 0)
		{
			direction = !negSourceCats.contains(conCat); // if negative conditional cat, direction is decrease	
			finalDiffWeight = conWeight;
		}

		// concat is void refcat is not   --> differential decrease (unless refcat is negative)
		else if (conWeight == 0)
		{
			direction = negSourceCats.contains(refCat); // if negative ref cat, direction is positive
			finalDiffWeight = refWeight;
		}

		// refcat is positive, concat is negative 	--> differential decrease by sum of weights
		else if (posSourceCats.contains(refCat) && negSourceCats.contains(conCat))
		{
			direction = false;
			finalDiffWeight = refWeight + conWeight;
		}

		// refcat is negative, concat is positive --> differential increase by sum of weights
		else if (negSourceCats.contains(refCat) && posSourceCats.contains(conCat))
		{
			direction = true;
			finalDiffWeight = refWeight + conWeight;
		}

		else
		{
			finalDiffWeight = refWeight - conWeight;
			if (finalDiffWeight < 0) // conWeight is higher than refWeight: differential increase
			{
				finalDiffWeight *= -1;
				direction = true;
			}
			else if (finalDiffWeight > 0) // refWeight is higher than conWeight: differential decrease
			{
				direction = false;
			}
			boolean refNeutral = !(negSourceCats.contains(refCat) || posSourceCats.contains(refCat));
			boolean conNeutral = !(negSourceCats.contains(conCat) || posSourceCats.contains(conCat));

			if ((refNeutral && !conNeutral) || (conNeutral && !refNeutral))
			{
				unspecified = true;
			}
		}
		if (finalDiffWeight <= weightCutoff)
		{
			return eg.getVoidEdge(final_symm);
		}

		// determine final type by added prefixes such as 'unspecified' and 'increase'
		String type = "";
		if (direction)
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
		else
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
		if (unspecified)
		{
			type += teo.getUnspecifiedPrefix();
		}
		type += baseType;
		diff_edge.setType(type);
		diff_edge.setWeight(finalDiffWeight);
		return diff_edge;
	}

	/**
	 * Method that defines the consensus edge from the corresponding edge categories in the reference and condition-specific networks.
	 * Returns an empty set when the edge should be deleted (i.e. not present in the consensus network).
	 * 
	 * An important parameter is supportingCutoff, which determines the amount of support needed for an edge to be included in the consensus network. If it equals the number of input networks, all networks need to agree on an edge.
	 * However, if it is smaller, e.g. 3 out of 4, there can be one 'outlier' network (potentially a different one for each calculated edge), allowing some noise in the input and creating more robust consensus networks.
	 * The supportingCutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting edges (networks).
	 * 
	 * @param edges the original edge definitions (should not be null or empty!)
	 * @param supports the network IDs which support the corresponding edges in the other list 
	 * @param supportingCutoff the minimal number of networks that need to support an edge in the consensus network
	 * @param weightCutoff the minimal value of a resulting edge for it to be included in the consensus network
	 * @param minOperator whether or not to take the minimum of the edge weights - if false, the maximum is taken
	 * 
	 * @return the edge definitions in the consensus network, or an empty set, but never null
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public Map<EdgeDefinition, Set<Integer>> getConsensusEdge(List<EdgeDefinition> edges, List<Set<Integer>> supports, int supportingCutoff, double weightCutoff, boolean minOperator) throws IllegalArgumentException
	{
		Map<EdgeDefinition, Set<Integer>> consensusEdges = new HashMap<EdgeDefinition, Set<Integer>>();

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

		// 2. GO THROUGH THE WHOLE ONTOLOGY TREE AND COLLECT THE RESULTS, USING THE SUPPORT CUTOFF PARAMETER //

		for (String cat : teo.getAllSourceCategories(true))
		{
			IntermediateComparison aff_result = affirmative_results.get(cat);
			boolean aff = false;
			if (aff_result != null && aff_result.getTotalSupport() >= supportingCutoff)
			{
				aff = true;
			}

			IntermediateComparison neg_result = negated_results.get(cat);
			boolean neg = false;
			if (neg_result != null && neg_result.getTotalSupport() >= supportingCutoff)
			{
				neg = true;
			}

			if (aff)
			{
				Map<EdgeDefinition, Set<Integer>> map = createAllEdges(aff_result, final_symm, false, supportingCutoff, weightCutoff, Double.POSITIVE_INFINITY, minOperator);

				for (EdgeDefinition consensusEdge : map.keySet())
				{
					Set<Integer> theseSupports = map.get(consensusEdge);
					consensusEdges.put(consensusEdge, theseSupports);
				}
			}
			if (neg)
			{
				Map<EdgeDefinition, Set<Integer>> map = createAllEdges(neg_result, final_symm, true, supportingCutoff, weightCutoff, Double.POSITIVE_INFINITY, minOperator);
				for (EdgeDefinition consensusEdge : map.keySet())
				{
					Set<Integer> theseSupports = map.get(consensusEdge);
					consensusEdges.put(consensusEdge, theseSupports);
				}
			}
		}

		return consensusEdges;
	}

}

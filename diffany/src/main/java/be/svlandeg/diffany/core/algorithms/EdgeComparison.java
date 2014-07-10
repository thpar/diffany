package be.svlandeg.diffany.core.algorithms;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

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
	 * @param eo the edge ontology defining the interaction semantics
	 */
	public EdgeComparison(TreeEdgeOntology teo)
	{
		this.teo = teo;
		eg = new EdgeGenerator();
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

		String firstParent = teo.retrieveFirstCommonParent(allEdges, true);
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
	 * Returns EdgeDefinition.getVoidEdge() when the edge should be deleted (i.e. not present in the overlapping network).
	 * 
	 * @param edges the original edge definitions (can contain EdgeDefinition.getVoidEdge()), should not be empty!
	 * @param noNetworks the number of original input networks
	 * @param no_cutoff the number of networks that need to have the overlap for it to be included
	 * @param weight_cutoff the minimal value of a resulting edge for it to be included in the overlapping network
	 * @param minOperator whether or not to take the minimum of the edge weights - if false, the maximum is taken

	 * @return the edge definition in the overlapping network, or EdgeDefinition.getVoidEdge() when there should be no such edge (never null).
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public EdgeDefinition getOverlapEdge(Collection<EdgeDefinition> edges, int noNetworks, int no_cutoff, double weight_cutoff, boolean minOperator) throws IllegalArgumentException
	{
		// TODO use no_cutoff properly
		
		if (edges == null || edges.isEmpty())
		{
			String errormsg = "The set of edges should not be null or empty!";
			throw new IllegalArgumentException(errormsg);
		}
		EdgeDefinition overlap_edge = eg.getDefaultEdge();
		int countEdges = edges.size();
		
		// 1. DETERMINE NEGATION AND SYMMETRY //
		int countNegated = 0;
		int countSymmetrical = 0;

		double minWeight = Double.MAX_VALUE;
		double maxWeight = Double.MIN_VALUE;

		for (EdgeDefinition e : edges)
		{
			if (e.isNegated())
			{
				countNegated++;
			}
			
			if (e.isSymmetrical())
			{
				countSymmetrical++;
			}

			double weight = e.getWeight();
			if (weight < minWeight)
			{
				minWeight = weight;
			}
			if (weight > maxWeight)
			{
				maxWeight = weight;
			}
		}

		boolean symm = countSymmetrical == countEdges;
		overlap_edge.makeSymmetrical(symm);
		
		// If there are less input edges than the cutoff requires, there will not be an overlap edge, return eg.getVoidEdge(symm)
		if (countEdges < no_cutoff)
		{
			return eg.getVoidEdge(symm);
		}

		// some are negated, some are not -> no overlap
		if (countNegated != 0 && countNegated != countEdges)
		{
			return eg.getVoidEdge(symm);
		}
		boolean overlapNegated = countNegated == countEdges;
		overlap_edge.makeNegated(overlapNegated);

		// 2. DETERMINE WEIGHT //

		// the overlapping weight is the minimum between the two, or the maximum
		// if specified as such
		double overlapWeight = minWeight;
		if (!minOperator)
		{
			overlapWeight = maxWeight;
		}
		if (overlapWeight <= weight_cutoff)
		{
			return eg.getVoidEdge(symm);
		}
		overlap_edge.setWeight(overlapWeight);

		// 3. DEFINE TYPE BY INSPECTING CHILDREN AND PARENTS //

		String firstCommonParent = teo.retrieveFirstCommonParent(edges, false);
		if (firstCommonParent == null)
		{
			// no category covers all of the edges
			return eg.getVoidEdge(symm);
		}

		if (!overlapNegated && firstCommonParent != null) //  the shared edge is the (first) common super class 
		{
			overlap_edge.setType(firstCommonParent);
			return overlap_edge;
		}

		String firstCommonChild = teo.retrieveFirstCommonChild(edges);
		if (firstCommonParent == null)
		{
			// no category covers all of the edges
			return eg.getVoidEdge(symm);
		}

		// the shared edge is the negation of the (first) common subclass, if there is one such
		overlap_edge.setType(firstCommonChild);
		return overlap_edge;
	}

}

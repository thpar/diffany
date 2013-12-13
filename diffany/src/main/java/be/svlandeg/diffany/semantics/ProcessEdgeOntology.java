package be.svlandeg.diffany.semantics;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This edge ontology deals with process edges such as ppi and ptm, their
 * interrelationships and their corresponding weights.
 * 
 * @author Sofie Van Landeghem
 */
public class ProcessEdgeOntology extends EdgeOntology
{

	private String negPrefix;
	private String posPrefix;

	/**
	 * Create a new ontology, defining the set of categories. and inserting
	 * After the constructor is called, default edge-category mappings should be
	 * inserted using addCategoryMapping!
	 */
	public ProcessEdgeOntology(String posPrefix, String negPrefix, Set<String> sourceCats)
	{
		super();
		this.negPrefix = negPrefix;
		this.posPrefix = posPrefix;
		addSourceCategories(sourceCats);
	}

	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, Set<EdgeDefinition> conEdges, double cutoff) throws IllegalArgumentException
	{
		EdgeDefinition diff_edge = new EdgeDefinition();
		Set<EdgeDefinition> conEdges2 = new HashSet<EdgeDefinition>();
		
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
		
		boolean refNeg = refEdge.isNegated();
		if (refNeg)
			refEdge = EdgeDefinition.getVoidEdge();
		
		diff_edge.makeNegated(false);	// a differential edge is never negated 
		
		boolean refSymm = refEdge.isSymmetrical();
		boolean diffSymm = refSymm && conSymm;
		diff_edge.makeSymmetrical(diffSymm);
		
		//////////// DETERMINE TYPE AND WEIGHT ////////////////
		
		String refCat = getSourceCategory(refEdge.getType());
		
		int countUp = 0;
		int countDown = 0;
		double minDiffWeight = Double.MAX_VALUE;
		double minCumulWeight = Double.MAX_VALUE;
		
		for (EdgeDefinition conEdge : conEdges2)
		{
			String conCat = getSourceCategory(conEdge.getType());
			double diffWeight = 0;
			double cumWeight = conEdge.getWeight() + refEdge.getWeight();
			
			if (refCat.equals(conCat)) 
			{
				diffWeight = conEdge.getWeight() - refEdge.getWeight();
				if (diffWeight < 0) // decrease
				{
					diffWeight *= -1;
					countDown++;
				}
				else if (diffWeight > 0)	// increase
				{
					countUp++;
				}
			}
			else if (refCat.equals(VOID_TYPE) && ! conCat.equals(VOID_TYPE))
			{
				countUp++;
				diffWeight = conEdge.getWeight();
			}

			else if (conCat.equals(VOID_TYPE) && ! refCat.equals(VOID_TYPE))
			{
				countDown++;
				diffWeight = refEdge.getWeight();
			}
			minDiffWeight = Math.min(minDiffWeight, diffWeight);
			minCumulWeight = Math.min(minCumulWeight, cumWeight);
		}
		// some are up, some are down -> no general differential edge
		if (countUp > 0 && countDown > 0)
		{
			return EdgeDefinition.getVoidEdge();
		}
		
		Boolean up = null;
		if (countUp == 0 && countDown == conEdges2.size())
		{
			// all differential edges are down
			up = false;
		}
		if (countDown == 0 && countUp == conEdges2.size())
		{
			// all differential edges are up
			up = true;
		}
		
		Map<String, Integer> conParents = retrieveFirstCommonParents(conEdges2);
		if (conParents.isEmpty()) 	
		{
			return EdgeDefinition.getVoidEdge();
		}
		String conParent = conParents.keySet().iterator().next();
		
		String type = "";
		double diffWeight = 0.0;
		if (up == null)
		{
			// fully changed category: take common parent of conEdges, if there is one
			// however, only do this when the two categories are not directly related (unspecified)
			
			int related1 = isSourceChildOf(refCat, conParent);
			int related2 = isSourceChildOf(conParent, refCat);
			
			if (related1 >= 0 || related2 >= 0)
			{
				return EdgeDefinition.getVoidEdge();
			}
			type = refCat + "_to_" + conParent;
			diffWeight = minCumulWeight;
		}
		else 
		{
			if (up)
			{
				type = posPrefix;
			} 
			else
			{
				type = negPrefix;
			}
			String baseType = refCat;
			if (refCat.equals(VOID_TYPE))
			{
				baseType = conParent;
			}
			type += baseType;
			diffWeight = minDiffWeight;
		}
		diff_edge.setType(type);
		
		if (diffWeight <= cutoff)
		{
			return EdgeDefinition.getVoidEdge();
		}
		diff_edge.setWeight(diffWeight);
		
		return diff_edge;
	}

}

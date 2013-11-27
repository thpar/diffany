package be.svlandeg.diffany.semantics;

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
	public ProcessEdgeOntology(String posPrefix, String negPrefix, Set<String> cats)
	{
		super();
		this.negPrefix = negPrefix;
		this.posPrefix = posPrefix;
		addCategories(cats);
	}

	@Override
	public EdgeDefinition getSharedEdge(EdgeDefinition refEdge, EdgeDefinition conEdge, double cutoff) throws IllegalArgumentException
	{
		EdgeDefinition shared_edge = new EdgeDefinition();

		String refCat = getCategory(refEdge.getType());
		String conCat = getCategory(conEdge.getType());

		if (refCat.equals(conCat))
		{
			// the shared weight is the minimum between the two
			double sharedWeight = Math.min(refEdge.getWeight(), conEdge.getWeight());

			if (sharedWeight <= cutoff)
			{
				return EdgeDefinition.getVoidEdge();
			}

			// the shared edge is only symmetrical if both original edges are
			boolean refSymm = refEdge.isSymmetrical();
			boolean conSymm = conEdge.isSymmetrical();
			boolean sharedSymm = refSymm && conSymm;

			// the shared edge is only negated if both original edges are
			boolean refNeg = refEdge.isNegated();
			boolean conNeg = conEdge.isNegated();
			boolean sharedNeg = refNeg && conNeg;

			shared_edge.setType(refEdge.getType());
			shared_edge.setWeight(sharedWeight);
			shared_edge.makeSymmetrical(sharedSymm);
			shared_edge.makeNegated(sharedNeg);
			return shared_edge;
		}
		return EdgeDefinition.getVoidEdge();
	}

	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, EdgeDefinition conEdge, double cutoff) throws IllegalArgumentException
	{
		EdgeDefinition diff_edge = new EdgeDefinition();

		String refCat = getCategory(refEdge.getType());
		String conCat = getCategory(conEdge.getType());

		boolean equalCats = refCat.equals(conCat);
		Boolean up = null;
		double diffWeight = 0;

		if (equalCats)
		{
			diffWeight = conEdge.getWeight() - refEdge.getWeight();
			up = true; // assume an increase
			if (diffWeight < 0) // it was a decrease!
			{
				diffWeight *= -1;
				up = false;
			}
		}

		if (refCat.equals(VOID_TYPE))
		{
			up = true;
			diffWeight = conEdge.getWeight();
		}

		if (conCat.equals(VOID_TYPE))
		{
			up = false;
			diffWeight = refEdge.getWeight();
		}

		if (up == null || diffWeight <= cutoff)
		{
			return EdgeDefinition.getVoidEdge();
		}

		String type = "";
		if (up)
		{
			type = posPrefix;
		} else
		{
			type = negPrefix;
		}
		type += conCat;

		diff_edge.setType(type);

		// the differential edge is only symmetrical if both original edges are
		// TODO is this correct, also for void types?
		boolean refSymm = refEdge.isSymmetrical();
		boolean conSymm = conEdge.isSymmetrical();
		boolean diffSymm = refSymm && conSymm;

		// the differential edge is only negated if both original edges are
		// TODO is this correct, also for void types?
		boolean refNeg = refEdge.isNegated();
		boolean conNeg = conEdge.isNegated();
		boolean diffNeg = refNeg && conNeg;

		diff_edge.setWeight(diffWeight);
		diff_edge.makeSymmetrical(diffSymm);
		diff_edge.makeNegated(diffNeg);
		return diff_edge;
	}

}

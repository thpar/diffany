package be.svlandeg.diffany.semantics;

import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * The up/down edge ontology can deal with up/down regulation and according differences in weight.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class UpDownEdgeOntology extends EdgeOntology
{

	protected String neg_diff_cat;
	protected String pos_diff_cat;

	public Set<String> posCats;
	public Set<String> negCats;
	public Set<String> neutralCats;

	/**
	 * Create a new ontology, defining pos/neg/neutral categories and inserting default edge-category mappings.
	 */
	public UpDownEdgeOntology(String pos_diff_cat, String neg_diff_cat)
	{
		this.neg_diff_cat = neg_diff_cat;
		this.pos_diff_cat = pos_diff_cat;
		removeAllCategoriesAndMappings();
		definePosCategories();
		defineNegCategories();
		defineNeutralCategories();
		insertDefaultMappings();
	}

	/**
	 * Define all the positive categories that are defined in this ontology.
	 */
	protected abstract void definePosCategories();

	/**
	 * Define all the negative categories that are defined in this ontology.
	 */
	protected abstract void defineNegCategories();

	/**
	 * Define all the neutral categories that are defined in this ontology.
	 */
	protected abstract void defineNeutralCategories();

	/**
	 * Provide a default mapping edge type to category mapping. 
	 */
	protected abstract void insertDefaultMappings();

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

			if (sharedWeight < cutoff)
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

		if (posCats.contains(refCat) && posCats.contains(conCat))
		{
			equalCats = true;
			up = true;
		}

		if (negCats.contains(refCat) && negCats.contains(conCat))
		{
			equalCats = true;
			up = false;
		}

		if (posCats.contains(refCat) && negCats.contains(conCat))
		{
			up = false;
			equalCats = false;
		}
		if (posCats.contains(refCat) && conCat.equals(VOID_TYPE))
		{
			up = false;
			equalCats = false;
		}
		if (refCat.equals(VOID_TYPE) && negCats.contains(conCat))
		{
			up = false;
			equalCats = false;
		}

		if (negCats.contains(refCat) && posCats.contains(conCat))
		{
			up = true;
			equalCats = false;
		}
		if (negCats.contains(refCat) && conCat.equals(VOID_TYPE))
		{
			up = true;
			equalCats = false;
		}
		if (refCat.equals(VOID_TYPE) && posCats.contains(conCat))
		{
			up = true;
			equalCats = false;
		}

		double diffWeight = 0;
		if (equalCats)
		{
			diffWeight = conEdge.getWeight() - refEdge.getWeight();

			// the weight has decreased within equal categories, which means the up/down direction changes
			if (diffWeight < 0)
			{
				diffWeight *= -1;
				up = !up;
			}
		}
		else
		{
			diffWeight = conEdge.getWeight() + refEdge.getWeight();
		}

		if (up == null || diffWeight < cutoff)
		{
			return EdgeDefinition.getVoidEdge();
		}

		if (up)
			diff_edge.setType(pos_diff_cat);
		if (!up)
			diff_edge.setType(neg_diff_cat);

		// the differential edge is only symmetrical if both original edges are
		// TODO is this correct?
		boolean refSymm = refEdge.isSymmetrical();
		boolean conSymm = conEdge.isSymmetrical();
		boolean diffSymm = refSymm && conSymm;

		// the differential edge is only negated if both original edges are
		// TODO is this correct?
		boolean refNeg = refEdge.isNegated();
		boolean conNeg = conEdge.isNegated();
		boolean diffNeg = refNeg && conNeg;

		diff_edge.setWeight(diffWeight);
		diff_edge.makeSymmetrical(diffSymm);
		diff_edge.makeNegated(diffNeg);
		return diff_edge;

	}
}

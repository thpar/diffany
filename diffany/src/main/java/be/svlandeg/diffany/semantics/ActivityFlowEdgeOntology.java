package be.svlandeg.diffany.semantics;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This edge ontology deals with general flow activity as up/down regulation and their corresponding weights.
 * TODO: deal with neutral/undefined flow?
 * 
 * @author Sofie Van Landeghem
 */
public class ActivityFlowEdgeOntology extends EdgeOntology
{

	protected String neg_diff_cat;
	protected String pos_diff_cat;

	public Set<String> posCats;
	public Set<String> negCats;

	/**
	 * Create a new ontology, defining pos/neg categories. and inserting
	 * After the constructor is called, default edge-category mappings should be inserted using addCategoryMapping!
	 */
	public ActivityFlowEdgeOntology(String pos_diff_cat, String neg_diff_cat, Set<String> other_pos_cats, Set<String> other_neg_cats)
	{
		super();
		this.neg_diff_cat = neg_diff_cat;
		this.pos_diff_cat = pos_diff_cat;
		definePosCategories(other_pos_cats);
		defineNegCategories(other_neg_cats);
	}

	/**
	 * Define all the positive categories that are defined in this ontology.
	 */
	protected void definePosCategories(Set<String> other_pos_cats)
	{
		posCats = new HashSet<String>();
		posCats.add(pos_diff_cat);

		for (String p : other_pos_cats)
		{
			posCats.add(p);
		}
		addCategories(posCats);
	}

	/**
	 * Define all the negative categories that are defined in this ontology.
	 */
	protected void defineNegCategories(Set<String> other_neg_cats)
	{
		negCats = new HashSet<String>();
		negCats.add(neg_diff_cat);

		for (String p : other_neg_cats)
		{
			negCats.add(p);
		}
		addCategories(negCats);
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

			// the weight has decreased within equal categories, which means the
			// up/down direction changes
			if (diffWeight < 0)
			{
				diffWeight *= -1;
				up = !up;
			}
		} else
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

package be.svlandeg.diffany.semantics;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This edge ontology deals with general flow activity as up/down regulation and their corresponding weights.
 * 
 * @author Sofie Van Landeghem
 */
public class ActivityFlowEdgeOntology extends EdgeOntology
{

	protected String neg_diff_cat;
	protected String pos_diff_cat;

	public Set<String> posCats;
	public Set<String> negCats;
	public Set<String> neutralCats;

	/**
	 * Create a new ontology, defining pos/neg categories. and inserting
	 * After the constructor is called, default edge-category mappings should be inserted using addCategoryMapping!
	 * 
	 * @param pos_diff_cat the string representing an increase in a differential network
	 * @param neg_diff_cat the string representing an decrease in a differential network
	 * @param other_pos_cats positive categories (e.g. upregulation)
	 * @param other_neg_cats negative categories (e.g. inhibition)
	 * @param neutral_cats neutral categories (e.g. regulation)
	 */
	public ActivityFlowEdgeOntology(String pos_diff_cat, String neg_diff_cat, Set<String> other_pos_cats, Set<String> other_neg_cats, Set<String> other_neutral_cats)
	{
		super();
		this.neg_diff_cat = neg_diff_cat;
		this.pos_diff_cat = pos_diff_cat;
		definePosCategories(other_pos_cats);
		defineNegCategories(other_neg_cats);
		defineNeutralCategories(other_neutral_cats);
	}

	/**
	 * Define all the positive categories that are defined in this ontology.
	 * @param other_pos_cats the positive categories
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
	 * @param other_neg_cats the negative categories
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
	
	/**
	 * Define all the neutral categories that are defined in this ontology.
	 * @param neutral_cats the neutral categories
	 */
	protected void defineNeutralCategories(Set<String> all_neutral_cats)
	{
		neutralCats = new HashSet<String>();

		for (String p : all_neutral_cats)
		{
			neutralCats.add(p);
		}
		addCategories(neutralCats);
	}


	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, EdgeDefinition conEdge, double cutoff) throws IllegalArgumentException
	{
		EdgeDefinition diff_edge = new EdgeDefinition();

		boolean refNeg = refEdge.isNegated();
		boolean conNeg = conEdge.isNegated();
		
		if (refNeg)
			refEdge = EdgeDefinition.getVoidEdge();
		
		if (conNeg)
			conEdge = EdgeDefinition.getVoidEdge();
			
		String refCat = getCategory(refEdge.getType());
		String conCat = getCategory(conEdge.getType());
		
		// the differential edge can not be calculated if either of the two is a neutral edge
		if (neutralCats.contains(refCat) || neutralCats.contains(conCat))
		{
			return EdgeDefinition.getVoidEdge();
		}
		

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
		if (equalCats)		// within the same categories, edge weights are substracted from eachother
		{
			diffWeight = conEdge.getWeight() - refEdge.getWeight();

			// the weight has decreased within equal categories, which means the
			// up/down direction changes
			if (diffWeight < 0)
			{
				diffWeight *= -1;
				up = !up;
			}
		} else	// within opposite categories, edge weights are summed up to get the differential weight
		{
			diffWeight = conEdge.getWeight() + refEdge.getWeight();
		}

		if (up == null || diffWeight <= cutoff)
		{
			return EdgeDefinition.getVoidEdge();
		}

		if (up)
			diff_edge.setType(pos_diff_cat);
		if (!up)
			diff_edge.setType(neg_diff_cat);
		
		boolean refSymm = refEdge.isSymmetrical();
		boolean conSymm = conEdge.isSymmetrical();
		boolean diffSymm = refSymm && conSymm;

		diff_edge.setWeight(diffWeight);
		diff_edge.makeSymmetrical(diffSymm);
		diff_edge.makeNegated(false);	// a differential edge is never negated 
		
		return diff_edge;
	}
}

package be.svlandeg.diffany.semantics;

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

	public Set<String> source_pos_cats;
	public Set<String> source_neg_cats;
	public Set<String> source_neutral_cats;

	/**
	 * Create a new ontology, defining pos/neg categories. and inserting
	 * After the constructor is called, default edge-category mappings should be inserted using addCategoryMapping!
	 * 
	 * @param pos_diff_cat the string representing an increase in a differential network
	 * @param neg_diff_cat the string representing an decrease in a differential network
	 * @param source_pos_cats positive categories in the input networks (e.g. upregulation)
	 * @param source_neg_cats negative categories in the input networks (e.g. inhibition)
	 * @param source_neutral_cats neutral categories in the input networks (e.g. regulation)
	 */
	public ActivityFlowEdgeOntology(String pos_diff_cat, String neg_diff_cat, Set<String> source_pos_cats, Set<String> source_neg_cats, Set<String> source_neutral_cats)
	{
		super();
		
		this.pos_diff_cat = pos_diff_cat;
		this.neg_diff_cat = neg_diff_cat;
		
		addDiffCategory(pos_diff_cat);
		addDiffCategory(neg_diff_cat);
		
		this.source_pos_cats = source_pos_cats;
		this.source_neg_cats = source_neg_cats;
		this.source_neutral_cats = source_neutral_cats;
		
		addSourceCategories(source_pos_cats);
		addSourceCategories(source_neg_cats);
		addSourceCategories(source_neutral_cats);
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
			
		String refCat = getSourceCategory(refEdge.getType());
		String conCat = getSourceCategory(conEdge.getType());
		
		// the differential edge can not be calculated if either of the two is a neutral edge
		if (source_neutral_cats.contains(refCat) || source_neutral_cats.contains(conCat))
		{
			return EdgeDefinition.getVoidEdge();
		}
		

		boolean equalCats = refCat.equals(conCat);

		Boolean up = null;

		if (source_pos_cats.contains(refCat) && source_pos_cats.contains(conCat))
		{
			equalCats = true;
			up = true;
		}

		if (source_neg_cats.contains(refCat) && source_neg_cats.contains(conCat))
		{
			equalCats = true;
			up = false;
		}

		if (source_pos_cats.contains(refCat) && source_neg_cats.contains(conCat))
		{
			up = false;
			equalCats = false;
		}
		if (source_pos_cats.contains(refCat) && conCat.equals(VOID_TYPE))
		{
			up = false;
			equalCats = false;
		}
		if (refCat.equals(VOID_TYPE) && source_neg_cats.contains(conCat))
		{
			up = false;
			equalCats = false;
		}

		if (source_neg_cats.contains(refCat) && source_pos_cats.contains(conCat))
		{
			up = true;
			equalCats = false;
		}
		if (source_neg_cats.contains(refCat) && conCat.equals(VOID_TYPE))
		{
			up = true;
			equalCats = false;
		}
		if (refCat.equals(VOID_TYPE) && source_pos_cats.contains(conCat))
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

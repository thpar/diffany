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
	public EdgeDefinition getSharedEdge(EdgeDefinition referenceEdge, EdgeDefinition conditionEdge) throws IllegalArgumentException
	{
		EdgeDefinition shared_edge = new EdgeDefinition();

		String refCat = getCategory(referenceEdge.getType());
		String conCat = getCategory(conditionEdge.getType());

		if (refCat.equals(conCat))
		{
			shared_edge.setType(referenceEdge.getType());
			shared_edge.setWeight(referenceEdge.getWeight());
			shared_edge.makeSymmetrical(referenceEdge.isSymmetrical());
			shared_edge.makeNegated(referenceEdge.isNegated());
			return shared_edge;
		}
		return EdgeDefinition.getVoidEdge();
	}

	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition referenceEdge, EdgeDefinition conditionEdge) throws IllegalArgumentException
	{
		EdgeDefinition diff_edge = new EdgeDefinition();

		String refCat = getCategory(referenceEdge.getType());
		String conCat = getCategory(conditionEdge.getType());

		Boolean up = null;

		if (posCats.contains(refCat) && negCats.contains(conCat))
		{
			up = false;
		}
		if (posCats.contains(refCat) && conCat.equals(VOID_TYPE))
		{
			up = false;
		}
		if (refCat.equals(VOID_TYPE) && negCats.contains(conCat))
		{
			up = false;
		}

		if (negCats.contains(refCat) && posCats.contains(conCat))
		{
			up = true;
		}
		if (negCats.contains(refCat) && conCat.equals(VOID_TYPE))
		{
			up = true;
		}
		if (refCat.equals(VOID_TYPE) && posCats.contains(conCat))
		{
			up = true;
		}

		if (up != null)
		{
			if (up)
				diff_edge.setType(pos_diff_cat);
			if (!up)
				diff_edge.setType(neg_diff_cat);
			diff_edge.setWeight(referenceEdge.getWeight());
			diff_edge.makeSymmetrical(referenceEdge.isSymmetrical());
			diff_edge.makeNegated(referenceEdge.isNegated());
			return diff_edge;
		}

		return EdgeDefinition.getVoidEdge();
	}
}

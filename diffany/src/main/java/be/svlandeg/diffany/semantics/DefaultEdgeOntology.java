package be.svlandeg.diffany.semantics;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This class creates a default ActivityFlowEdgeOntology that can be used to
 * give initial suggestions to the user concerning edge type-to-category
 * mappings.
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultEdgeOntology extends EdgeOntology
{

	private ActivityFlowEdgeOntology afOntology;
	private ProcessEdgeOntology prOntology;

	/**
	 * Create a default edge ontology, with generic up/down edge categories and
	 * type-category mappings.
	 */
	public DefaultEdgeOntology()
	{
		String pos_cat = "increase";
		Set<String> other_pos_cats = definePosCategories();

		String neg_cat = "decrease";
		Set<String> other_neg_cats = defineNegCategories();

		Set<String> process_cats = defineProcessCategories();

		afOntology = new ActivityFlowEdgeOntology(pos_cat, neg_cat, other_pos_cats, other_neg_cats);
		prOntology = new ProcessEdgeOntology(pos_cat + "_", neg_cat + "_", process_cats);

		insertDefaultMappings();

	}

	/**
	 * Define the process categories to be used in the ProcessEdgeOntology
	 * 
	 * @return the defined process categories
	 */
	protected Set<String> defineProcessCategories()
	{
		Set<String> cats = new HashSet<String>();
		cats.add("ppi");
		return cats;
	}

	/**
	 * Define the positive categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined positive categories
	 */
	protected Set<String> definePosCategories()
	{
		Set<String> other_pos_cats = new HashSet<String>();
		other_pos_cats.add("pos_regulation");
		return other_pos_cats;
	}

	/**
	 * Define the negative categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined negative categories
	 */
	protected Set<String> defineNegCategories()
	{
		Set<String> other_neg_cats = new HashSet<String>();
		other_neg_cats.add("neg_regulation");
		return other_neg_cats;
	}

	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, EdgeDefinition conEdge, double cutoff) throws IllegalArgumentException
	{
		String refType = refEdge.getType();
		String conType = conEdge.getType();
		if (afOntology.isDefined(refType) && afOntology.isDefined(conType))
		{
			return afOntology.getDifferentialEdge(refEdge, conEdge, cutoff);
		}
		if (prOntology.isDefined(refType) && prOntology.isDefined(conType))
		{
			return prOntology.getDifferentialEdge(refEdge, conEdge, cutoff);
		}
		String errormsg = "The two types '" + refType + "' and '" + conType + "' can not be compared in this edge ontology!";
		throw new IllegalArgumentException(errormsg);
	}

	@Override
	public EdgeDefinition getSharedEdge(EdgeDefinition refEdge, EdgeDefinition conEdge, double cutoff) throws IllegalArgumentException
	{
		String refType = refEdge.getType();
		String conType = conEdge.getType();
		if (afOntology.isDefined(refType) && afOntology.isDefined(conType))
		{
			return afOntology.getSharedEdge(refEdge, conEdge, cutoff);
		}
		if (prOntology.isDefined(refType) && prOntology.isDefined(conType))
		{
			return prOntology.getSharedEdge(refEdge, conEdge, cutoff);
		}
		String errormsg = "The two types '" + refType + "' and '" + conType + "' can not be compared in this edge ontology!";
		throw new IllegalArgumentException(errormsg);
	}

	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void insertDefaultMappings()
	{
		boolean overwrite = false;

		// PPI category and common synonyms

		prOntology.addCategoryMapping("ppi", "ppi", overwrite);
		prOntology.addCategoryMapping("protein-protein interaction", "ppi", overwrite);
		prOntology.addCategoryMapping("proteinprotein interaction", "ppi", overwrite);
		prOntology.addCategoryMapping("proteinproteininteraction", "ppi", overwrite);
		prOntology.addCategoryMapping("binds", "ppi", overwrite);
		prOntology.addCategoryMapping("bind", "ppi", overwrite);
		prOntology.addCategoryMapping("binding", "ppi", overwrite);

		// regulation category and common synonyms
		// addCategoryMapping("regulation", "regulation", overwrite);
		// addCategoryMapping("regulates", "regulation", overwrite);
		// addCategoryMapping("regulate", "regulation", overwrite);

		// positive regulation category and common synonyms
		afOntology.addCategoryMapping("positive regulation", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positive_regulation", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positive-regulation", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positiveregulation", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positive reg", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positive-reg", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positive_reg", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positivereg", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("pos regulation", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("pos_regulation", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("pos-regulation", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("posregulation", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("pos reg", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("pos-reg", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("pos_reg", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("posreg", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("pos", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positive", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positively regulates", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positively_regulates", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positively-regulates", "pos_regulation", overwrite);
		afOntology.addCategoryMapping("positivelyregulates", "pos_regulation", overwrite);

		// positive regulation category and common synonyms
		afOntology.addCategoryMapping("negative regulation", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negative_regulation", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negative-regulation", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negativeregulation", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negative reg", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negative_reg", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negative-reg", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negativereg", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("neg regulation", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("neg_regulation", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("neg-regulation", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negregulation", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("neg reg", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("neg_reg", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("neg-reg", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negreg", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("neg", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negative", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negatively regulates", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negatively_regulates", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negatively-regulates", "neg_regulation", overwrite);
		afOntology.addCategoryMapping("negativelyregulates", "neg_regulation", overwrite);
	}

}

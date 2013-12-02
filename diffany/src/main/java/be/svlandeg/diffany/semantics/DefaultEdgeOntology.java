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
		
		Set<String> neutral_cats = defineNeutralCategories();
		
		afOntology = new ActivityFlowEdgeOntology(pos_cat, neg_cat, other_pos_cats, other_neg_cats, neutral_cats);

		Set<String> process_cats = defineProcessCategories();

		prOntology = new ProcessEdgeOntology(pos_cat + "_", neg_cat + "_", process_cats);

		insertDefaultParents();
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
		cats.add("ptm");
		cats.add("phosphorylation");
		cats.add("ubiquitination");
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
		other_pos_cats.add("positive_regulation");
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
		other_neg_cats.add("negative_regulation");
		return other_neg_cats;
	}
	
	/**
	 * Define the neutral categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined negative categories
	 */
	protected Set<String> defineNeutralCategories()
	{
		Set<String> neutral_cats = new HashSet<String>();
		neutral_cats.add("regulation");
		return neutral_cats;
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
	protected void insertDefaultParents()
	{
		afOntology.putParent("positive_regulation", "regulation");
		afOntology.putParent("negative_regulation", "regulation");
		
		prOntology.putParent("phosphorylation", "ptm");
		prOntology.putParent("ubiquitination", "ptm");
	}

	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void insertDefaultMappings()
	{
		boolean overwrite = false;

		// PPI category and common synonyms
		
		prOntology.addCategoryMapping("phosphorylation", "phosphorylation", overwrite);
		prOntology.addCategoryMapping("phosphorylates", "phosphorylation", overwrite);
		prOntology.addCategoryMapping("phosphorylate", "phosphorylation", overwrite);
		
		prOntology.addCategoryMapping("ubiquitination", "ubiquitination", overwrite);
		prOntology.addCategoryMapping("ubiquitinates", "ubiquitination", overwrite);
		prOntology.addCategoryMapping("ubiquitinate", "ubiquitination", overwrite);
		
		prOntology.addCategoryMapping("ptm", "ptm", overwrite);

		prOntology.addCategoryMapping("ppi", "ppi", overwrite);
		prOntology.addCategoryMapping("protein-protein interaction", "ppi", overwrite);
		prOntology.addCategoryMapping("proteinprotein interaction", "ppi", overwrite);
		prOntology.addCategoryMapping("proteinproteininteraction", "ppi", overwrite);
		prOntology.addCategoryMapping("binds", "ppi", overwrite);
		prOntology.addCategoryMapping("bind", "ppi", overwrite);
		prOntology.addCategoryMapping("binding", "ppi", overwrite);

		// regulation category and common synonyms
		afOntology.addCategoryMapping("regulation", "regulation", overwrite);
		afOntology.addCategoryMapping("regulates", "regulation", overwrite);
		afOntology.addCategoryMapping("regulate", "regulation", overwrite);
		afOntology.addCategoryMapping("influence", "regulation", overwrite);
		afOntology.addCategoryMapping("influences", "regulation", overwrite);
		afOntology.addCategoryMapping("effect", "regulation", overwrite);
		afOntology.addCategoryMapping("effects", "regulation", overwrite);

		// positive regulation category and common synonyms
		afOntology.addCategoryMapping("positive regulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positive_regulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positive-regulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positiveregulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positive reg", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positive-reg", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positive_reg", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positivereg", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("pos regulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("pos_regulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("pos-regulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("posregulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("pos reg", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("pos-reg", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("pos_reg", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("posreg", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("pos", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positive", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positively regulates", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positively_regulates", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positively-regulates", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("positivelyregulates", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("upregulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("up-regulation", "positive_regulation", overwrite);
		afOntology.addCategoryMapping("up regulation", "positive_regulation", overwrite);

		// positive regulation category and common synonyms
		afOntology.addCategoryMapping("negative regulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negative_regulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negative-regulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negativeregulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negative reg", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negative_reg", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negative-reg", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negativereg", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("neg regulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("neg_regulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("neg-regulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negregulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("neg reg", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("neg_reg", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("neg-reg", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negreg", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("neg", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negative", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negatively regulates", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negatively_regulates", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negatively-regulates", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("negativelyregulates", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("downregulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("down-regulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("down regulation", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("inhibits", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("inhibit", "negative_regulation", overwrite);
		afOntology.addCategoryMapping("inhibition", "negative_regulation", overwrite);
	}

}

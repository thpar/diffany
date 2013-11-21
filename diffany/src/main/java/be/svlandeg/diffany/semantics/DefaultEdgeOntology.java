package be.svlandeg.diffany.semantics;

import java.util.HashSet;

/**
 * This class creates a default UpDownEdgeOntology that can be used to give initial
 * suggestions to the user concerning edge type-to-category mappings.
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultEdgeOntology extends UpDownEdgeOntology
{

	/**
	 * Create a default edge ontology, with generic up/down edge categories and type-category mappings.
	 */
	public DefaultEdgeOntology()
	{
		super("increase", "decrease");
	}

	@Override
	protected void definePosCategories()
	{
		posCats = new HashSet<String>();

		posCats.add("pos_regulation");
		posCats.add(pos_diff_cat);

		addCategories(posCats);
	}

	@Override
	protected void defineNegCategories()
	{
		negCats = new HashSet<String>();

		negCats.add("neg_regulation");
		negCats.add(neg_diff_cat);

		addCategories(negCats);
	}

	@Override
	protected void defineNeutralCategories()
	{
		neutralCats = new HashSet<String>();
		neutralCats.add("regulation");
		neutralCats.add("ppi");

		addCategories(neutralCats);
	}

	@Override
	protected void insertDefaultMappings()
	{
		boolean overwrite = false;

		// PPI category and common synonyms
		addCategoryMapping("ppi", "ppi", overwrite);
		addCategoryMapping("protein-protein interaction", "ppi", overwrite);
		addCategoryMapping("proteinprotein interaction", "ppi", overwrite);
		addCategoryMapping("proteinproteininteraction", "ppi", overwrite);
		addCategoryMapping("binds", "ppi", overwrite);
		addCategoryMapping("bind", "ppi", overwrite);
		addCategoryMapping("binding", "ppi", overwrite);

		// regulation category and common synonyms
		addCategoryMapping("regulation", "regulation", overwrite);
		addCategoryMapping("regulates", "regulation", overwrite);
		addCategoryMapping("regulate", "regulation", overwrite);

		// positive regulation category and common synonyms
		addCategoryMapping("positive regulation", "pos_regulation", overwrite);
		addCategoryMapping("positive_regulation", "pos_regulation", overwrite);
		addCategoryMapping("positive-regulation", "pos_regulation", overwrite);
		addCategoryMapping("positiveregulation", "pos_regulation", overwrite);
		addCategoryMapping("positive reg", "pos_regulation", overwrite);
		addCategoryMapping("positive-reg", "pos_regulation", overwrite);
		addCategoryMapping("positive_reg", "pos_regulation", overwrite);
		addCategoryMapping("positivereg", "pos_regulation", overwrite);
		addCategoryMapping("pos regulation", "pos_regulation", overwrite);
		addCategoryMapping("pos_regulation", "pos_regulation", overwrite);
		addCategoryMapping("pos-regulation", "pos_regulation", overwrite);
		addCategoryMapping("posregulation", "pos_regulation", overwrite);
		addCategoryMapping("pos reg", "pos_regulation", overwrite);
		addCategoryMapping("pos-reg", "pos_regulation", overwrite);
		addCategoryMapping("pos_reg", "pos_regulation", overwrite);
		addCategoryMapping("posreg", "pos_regulation", overwrite);
		addCategoryMapping("pos", "pos_regulation", overwrite);
		addCategoryMapping("positive", "pos_regulation", overwrite);
		addCategoryMapping("positively regulates", "pos_regulation", overwrite);
		addCategoryMapping("positively_regulates", "pos_regulation", overwrite);
		addCategoryMapping("positively-regulates", "pos_regulation", overwrite);
		addCategoryMapping("positivelyregulates", "pos_regulation", overwrite);

		// positive regulation category and common synonyms
		addCategoryMapping("negative regulation", "neg_regulation", overwrite);
		addCategoryMapping("negative_regulation", "neg_regulation", overwrite);
		addCategoryMapping("negative-regulation", "neg_regulation", overwrite);
		addCategoryMapping("negativeregulation", "neg_regulation", overwrite);
		addCategoryMapping("negative reg", "neg_regulation", overwrite);
		addCategoryMapping("negative_reg", "neg_regulation", overwrite);
		addCategoryMapping("negative-reg", "neg_regulation", overwrite);
		addCategoryMapping("negativereg", "neg_regulation", overwrite);
		addCategoryMapping("neg regulation", "neg_regulation", overwrite);
		addCategoryMapping("neg_regulation", "neg_regulation", overwrite);
		addCategoryMapping("neg-regulation", "neg_regulation", overwrite);
		addCategoryMapping("negregulation", "neg_regulation", overwrite);
		addCategoryMapping("neg reg", "neg_regulation", overwrite);
		addCategoryMapping("neg_reg", "neg_regulation", overwrite);
		addCategoryMapping("neg-reg", "neg_regulation", overwrite);
		addCategoryMapping("negreg", "neg_regulation", overwrite);
		addCategoryMapping("neg", "neg_regulation", overwrite);
		addCategoryMapping("negative", "neg_regulation", overwrite);
		addCategoryMapping("negatively regulates", "neg_regulation", overwrite);
		addCategoryMapping("negatively_regulates", "neg_regulation", overwrite);
		addCategoryMapping("negatively-regulates", "neg_regulation", overwrite);
		addCategoryMapping("negativelyregulates", "neg_regulation", overwrite);
	}

}

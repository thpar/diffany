package be.svlandeg.diffany.semantics;


/**
 * This class creates a default TrippleEdgeOntology that can be used to give initial
 * suggestions to the user concerning edge type-to-category mappings.
 * 
 * @author Sofie Van Landeghem
 */
public class DefaultEdgeOntology extends TrippleEdgeOntology
{

	/**
	 * Create a default edge ontology, with generic edge categories and type-category mappings.
	 */
	public DefaultEdgeOntology()
	{
		super();
		removeAllCategoriesAndMappings();
		defineCategories();
		insertDefaultMappings();
		insertDefaultTripples();
	}

	/**
	 * Define all the categories that are defined in this ontology.
	 */
	protected void defineCategories()
	{
		allCategories.put("pos_regulation", false);
		allCategories.put("neg_regulation", false);
		allCategories.put("regulation", false);

		allCategories.put("decrease", false);
		allCategories.put("increase", false);

		allCategories.put("ppi", true);
	}

	/**
	 * Create default tripples, translating the most common edge categories to differential edges
	 * TODO revise current mappings from literature and possibly add more fancy stuff 
	 */
	private void insertDefaultTripples()
	{
		addTripple("pos_regulation", "neg_regulation", "decrease");
		addTripple("pos_regulation", "regulation", "decrease");
		addTripple("pos_regulation", VOID_EDGE, "decrease");
		addTripple(VOID_EDGE, "neg_regulation", "decrease");

		addTripple("regulation", "neg_regulation", "decrease");
		addTripple("regulation", VOID_EDGE, "decrease");
		
		addTripple("regulation", "pos_regulation", "increase");
		addTripple(VOID_EDGE, "regulation", "increase");

		addTripple("neg_regulation", "regulation", "increase");
		addTripple("neg_regulation", "pos_regulation", "increase");
		addTripple("neg_regulation", VOID_EDGE, "increase");
		addTripple(VOID_EDGE, "pos_regulation", "increase");
	}

	/**
	 * Default mapping. In TrippleEdgeOntology, upper/lower casing is not taken into
	 * account, so everything can be lower case.
	 * TODO: add more default mappings!
	 */
	private void insertDefaultMappings()
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

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
	 * Create a default edge ontology, with generic edge categories and type-category
	 * mappings.
	 */
	public DefaultEdgeOntology()
	{
		super();
		insertDefaultMappings();
		insertDefaultTripples();
	}

	
	/**
	 * Create default tripples, translating the most common edge categories to differential edges
	 * TODO revise current mappings from literature and possibly add more fancy stuff 
	 */
	public void insertDefaultTripples()
	{
		addTripple("pos_regulation","neg_regulation", "negative");
		addTripple("pos_regulation","regulation", "negative");
		
		addTripple("regulation","neg_regulation", "negative");
		addTripple("regulation","pos_regulation", "positive");
		
		addTripple("neg_regulation","regulation", "positive");
		addTripple("neg_regulation","pos_regulation", "positive");
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
		addCategoryMapping("ppi", "ppi", true, overwrite);
		addCategoryMapping("protein-protein interaction", "ppi", true,  overwrite);
		addCategoryMapping("proteinprotein interaction", "ppi", true,  overwrite);
		addCategoryMapping("proteinproteininteraction", "ppi", true,  overwrite);
		addCategoryMapping("binds", "ppi", true,  overwrite);
		addCategoryMapping("bind", "ppi", true,  overwrite);
		addCategoryMapping("binding", "ppi", true,  overwrite);

		// regulation category and common synonyms
		addCategoryMapping("regulation", "regulation", false, overwrite);
		addCategoryMapping("regulates", "regulation", false, overwrite);
		addCategoryMapping("regulate", "regulation", false, overwrite);

		// positive regulation category and common synonyms
		addCategoryMapping("positive regulation", "pos_regulation", false, overwrite);
		addCategoryMapping("positive_regulation", "pos_regulation", false, overwrite);
		addCategoryMapping("positive-regulation", "pos_regulation", false, overwrite);
		addCategoryMapping("positiveregulation", "pos_regulation", false, overwrite);
		addCategoryMapping("positive reg", "pos_regulation", false, overwrite);
		addCategoryMapping("positive-reg", "pos_regulation", false, overwrite);
		addCategoryMapping("positive_reg", "pos_regulation", false, overwrite);
		addCategoryMapping("positivereg", "pos_regulation", false, overwrite);
		addCategoryMapping("pos regulation", "pos_regulation", false, overwrite);
		addCategoryMapping("pos_regulation", "pos_regulation", false, overwrite);
		addCategoryMapping("pos-regulation", "pos_regulation", false, overwrite);
		addCategoryMapping("posregulation", "pos_regulation", false, overwrite);
		addCategoryMapping("pos reg", "pos_regulation", false, overwrite);
		addCategoryMapping("pos-reg", "pos_regulation", false, overwrite);
		addCategoryMapping("pos_reg", "pos_regulation", false, overwrite);
		addCategoryMapping("posreg", "pos_regulation", false, overwrite);
		addCategoryMapping("pos", "pos_regulation", false, overwrite);
		addCategoryMapping("positively regulates", "pos_regulation", false, overwrite);
		addCategoryMapping("positively_regulates", "pos_regulation", false, overwrite);
		addCategoryMapping("positively-regulates", "pos_regulation", false, overwrite);
		addCategoryMapping("positivelyregulates", "pos_regulation", false, overwrite);

		// positive regulation category and common synonyms
		addCategoryMapping("negative regulation", "neg_regulation", false, overwrite);
		addCategoryMapping("negative_regulation", "neg_regulation", false, overwrite);
		addCategoryMapping("negative-regulation", "neg_regulation", false, overwrite);
		addCategoryMapping("negativeregulation", "neg_regulation", false, overwrite);
		addCategoryMapping("negative reg", "neg_regulation", false, overwrite);
		addCategoryMapping("negative_reg", "neg_regulation", false, overwrite);
		addCategoryMapping("negative-reg", "neg_regulation", false, overwrite);
		addCategoryMapping("negativereg", "neg_regulation", false, overwrite);
		addCategoryMapping("neg regulation", "neg_regulation", false, overwrite);
		addCategoryMapping("neg_regulation", "neg_regulation", false, overwrite);
		addCategoryMapping("neg-regulation", "neg_regulation", false, overwrite);
		addCategoryMapping("negregulation", "neg_regulation", false, overwrite);
		addCategoryMapping("neg reg", "neg_regulation", false, overwrite);
		addCategoryMapping("neg_reg", "neg_regulation", false, overwrite);
		addCategoryMapping("neg-reg", "neg_regulation", false, overwrite);
		addCategoryMapping("negreg", "neg_regulation", false, overwrite);
		addCategoryMapping("neg", "neg_regulation", false, overwrite);
		addCategoryMapping("negatively regulates", "neg_regulation", false, overwrite);
		addCategoryMapping("negatively_regulates", "neg_regulation", false, overwrite);
		addCategoryMapping("negatively-regulates", "neg_regulation", false, overwrite);
		addCategoryMapping("negativelyregulates", "neg_regulation", false, overwrite);
	}

}

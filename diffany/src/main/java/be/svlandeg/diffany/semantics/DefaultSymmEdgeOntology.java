package be.svlandeg.diffany.semantics;

/**
 * This is a simplified version of the DefaultEdgeOntology, considering every interaction as symmetrical.
 *
 * @author Sofie Van Landeghem
 */
public class DefaultSymmEdgeOntology extends DefaultEdgeOntology
{
	
	/**
	 * Create a default edge ontology, with generic edge categories and type-category mappings.
	 */
	public DefaultSymmEdgeOntology()
	{
		super();
	}
	
	/**
	 * Define all the categories that are defined in this ontology.
	 */
	protected void defineCategories()
	{
		allCategories.put("pos_regulation", true);
		allCategories.put("neg_regulation", true);
		allCategories.put("regulation", true);

		allCategories.put("increase", true);
		allCategories.put("decrease", true);

		allCategories.put("ppi", true);
	}


}

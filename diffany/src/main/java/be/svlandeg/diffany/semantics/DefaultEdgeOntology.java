package be.svlandeg.diffany.semantics;

import java.awt.Paint;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This class holds both an ActivityFlowEdgeOntology and a ProcessEdgeOntology. 
 * It can be used to give initial suggestions to the user concerning edge type-to-category mappings.
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
		Set<String> pos_cats = definePosCategories();

		String neg_cat = "decrease";
		Set<String> neg_cats = defineNegCategories();
		
		Set<String> neutral_cats = defineNeutralCategories();
		
		afOntology = new ActivityFlowEdgeOntology(pos_cat, neg_cat, pos_cats, neg_cats, neutral_cats);

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
		cats.add("methylation");
		return cats;
	}

	/**
	 * Define the positive categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined positive categories
	 */
	protected Set<String> definePosCategories()
	{
		Set<String> pos_cats = new HashSet<String>();
		pos_cats.add("positive_regulation");
		return pos_cats;
	}

	/**
	 * Define the negative categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined negative categories
	 */
	protected Set<String> defineNegCategories()
	{
		Set<String> neg_cats = new HashSet<String>();
		neg_cats.add("negative_regulation");
		return neg_cats;
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
	
	/**
	 * Get the semantic category of a certain edge type. Matching is done independent of upper/lower casing.
	 * Both the internal ActivityFlowEdgeOntology as well as the ProcessEdgeOntology are consulted for the mapping.
	 * 
	 * @param edgeType the original type of the edge in a network
	 * @return the semantic category of that edge or null if it is not mapped in this ontology
	 */
	public String getSourceCategory(String edgeType)
	{
		String hit = afOntology.getSourceCategory(edgeType);
		if (hit == null)
		{
			hit = prOntology.getSourceCategory(edgeType);
		}
		return hit;
	}
	
	/**
	 * Return all categories present in this ontology, 
	 * consulting both the internal ActivityFlowEdgeOntology as well as the ProcessEdgeOntology.
	 * @return the set of categories mapped in this ontology.
	 */
	public Set<String> getAllSourceCategories()
	{
		Set<String> allCats = afOntology.getAllSourceCategories();
		allCats.addAll(prOntology.getAllSourceCategories());
		return allCats;
	}
	
	/**
	 * Return all categories present in this ontology, 
	 * consulting both the internal ActivityFlowEdgeOntology as well as the ProcessEdgeOntology.
	 * @return the set of categories mapped in this ontology.
	 */
	public Set<String> getAllDiffCategories()
	{
		Set<String> allCats = afOntology.getAllDiffCategories();
		allCats.addAll(prOntology.getAllDiffCategories());
		return allCats;
	}

	@Override
	public EdgeDefinition getDifferentialEdge(EdgeDefinition refEdge, Set<EdgeDefinition> conEdges, double cutoff) throws IllegalArgumentException
	{
		String refType = refEdge.getType();
		
		int conCount = conEdges.size();
		int afCount = 0;
		int prCount = 0;
		
		for (EdgeDefinition conEdge : conEdges)
		{
			String conType = conEdge.getType();
			if (afOntology.isDefinedSource(conType))
				afCount++;
			if (prOntology.isDefinedSource(conType))
				prCount++;
		}
		
		if (afOntology.isDefinedSource(refType) && afCount == conCount)
		{
			return afOntology.getDifferentialEdge(refEdge, conEdges, cutoff);
		}
		if (prOntology.isDefinedSource(refType) && prCount == conCount)
		{
			return prOntology.getDifferentialEdge(refEdge, conEdges, cutoff);
		}
		String errormsg = "The types '" + refType + "' and '" + conEdges.size() + "' conditional types could not be compared in this edge ontology!";
		throw new IllegalArgumentException(errormsg);
	}

	@Override
	public EdgeDefinition getOverlapEdge(Set<EdgeDefinition> edges, double cutoff, boolean minOperator) throws IllegalArgumentException
	{
		int conCount = edges.size();
		int afCount = 0;
		int prCount = 0;
		
		for (EdgeDefinition conEdge : edges)
		{
			String conType = conEdge.getType();
			if (afOntology.isDefinedSource(conType))
				afCount++;
			if (prOntology.isDefinedSource(conType))
				prCount++;
		}
		
		if (afCount == conCount)
		{
			return afOntology.getOverlapEdge(edges, cutoff, minOperator);
		}
		if (prCount == conCount)
		{
			return prOntology.getOverlapEdge(edges, cutoff, minOperator);
		}
		String errormsg = "The edge types could not be compared in this edge ontology!";
		throw new IllegalArgumentException(errormsg);
	}
	
	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void insertDefaultParents()
	{
		afOntology.putSourceParent("positive_regulation", "regulation");
		afOntology.putSourceParent("negative_regulation", "regulation");
		
		prOntology.putSourceParent("phosphorylation", "ptm");
		prOntology.putSourceParent("ubiquitination", "ptm");
		prOntology.putSourceParent("methylation", "ptm");
	}

	/**
	 * Provide a default mapping edge type to category mapping.
	 */
	protected void insertDefaultMappings()
	{
		boolean overwrite = false;

		// PPI category and common synonyms
		
		prOntology.addSourceCategoryMapping("methylation", "methylation", overwrite);
		prOntology.addSourceCategoryMapping("methylates", "methylation", overwrite);
		prOntology.addSourceCategoryMapping("methylate", "methylation", overwrite);
		
		prOntology.addSourceCategoryMapping("phosphorylation", "phosphorylation", overwrite);
		prOntology.addSourceCategoryMapping("phosphorylates", "phosphorylation", overwrite);
		prOntology.addSourceCategoryMapping("phosphorylate", "phosphorylation", overwrite);
		
		prOntology.addSourceCategoryMapping("ubiquitination", "ubiquitination", overwrite);
		prOntology.addSourceCategoryMapping("ubiquitinates", "ubiquitination", overwrite);
		prOntology.addSourceCategoryMapping("ubiquitinate", "ubiquitination", overwrite);
		
		prOntology.addSourceCategoryMapping("ptm", "ptm", overwrite);

		prOntology.addSourceCategoryMapping("ppi", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("protein-protein interaction", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("proteinprotein interaction", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("proteinproteininteraction", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("binds", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("bind", "ppi", overwrite);
		prOntology.addSourceCategoryMapping("binding", "ppi", overwrite);

		// regulation category and common synonyms
		afOntology.addSourceCategoryMapping("regulation", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("regulates", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("regulate", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("influence", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("influences", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("effect", "regulation", overwrite);
		afOntology.addSourceCategoryMapping("effects", "regulation", overwrite);

		// positive regulation category and common synonyms
		afOntology.addSourceCategoryMapping("positive regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive_regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive-regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positiveregulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive-reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive_reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positivereg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos_regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos-regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("posregulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos-reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos_reg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("posreg", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("pos", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positive", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positively regulates", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positively_regulates", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positively-regulates", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("positivelyregulates", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("upregulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("up-regulation", "positive_regulation", overwrite);
		afOntology.addSourceCategoryMapping("up regulation", "positive_regulation", overwrite);

		// negative regulation category and common synonyms
		afOntology.addSourceCategoryMapping("negative regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative_regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative-regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negativeregulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative_reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative-reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negativereg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg_regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg-regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negregulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg_reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg-reg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negreg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("neg", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negative", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negatively regulates", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negatively_regulates", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negatively-regulates", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("negativelyregulates", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("downregulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("down-regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("down regulation", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("inhibits", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("inhibit", "negative_regulation", overwrite);
		afOntology.addSourceCategoryMapping("inhibition", "negative_regulation", overwrite);
	}

	@Override
	public Paint getDifferentialEdgeStyle(EdgeOntology eo, String category)
	{
		if (afOntology.isDefinedCategory(category))
			return afOntology.getDifferentialEdgeStyle(eo, category);
		return prOntology.getDifferentialEdgeStyle(eo, category);
	}

	@Override
	public Paint getSourceEdgeStyle(EdgeOntology eo, String edgeType)
	{
		if (afOntology.isDefinedSource(edgeType))
			return afOntology.getSourceEdgeStyle(eo, edgeType);
		return prOntology.getSourceEdgeStyle(eo, edgeType);
	}

}

package be.svlandeg.diffany.semantics;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.EdgeDefinition;
import be.svlandeg.diffany.concepts.VisualEdgeStyle.ArrowHead;

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
		Map<String, Boolean> pos_cats = definePosCategories();

		String neg_cat = "decrease";
		Map<String, Boolean> neg_cats = defineNegCategories();
		
		Map<String, Boolean> neutral_cats = defineNeutralCategories();
		
		afOntology = new ActivityFlowEdgeOntology(pos_cat, neg_cat, pos_cats, neg_cats, neutral_cats);

		Map<String, Boolean> process_cats = defineProcessCategories();

		prOntology = new ProcessEdgeOntology(pos_cat + "_", neg_cat + "_", process_cats);

		insertDefaultParents();
		insertDefaultMappings();
		addDefaultPaintedParents();
	}

	/**
	 * Define the process categories to be used in the ProcessEdgeOntology
	 * 
	 * @return the defined process categories
	 */
	protected Map<String, Boolean> defineProcessCategories()
	{
		Map<String, Boolean> cats = new HashMap<String, Boolean>();
		cats.put("ppi", true);
		cats.put("ptm", false);
		cats.put("phosphorylation", false);
		cats.put("ubiquitination", false);
		cats.put("methylation", false);
		return cats;
	}

	/**
	 * Define the positive categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined positive categories
	 */
	protected Map<String, Boolean> definePosCategories()
	{
		Map<String, Boolean> pos_cats = new HashMap<String, Boolean>();
		pos_cats.put("positive_regulation", false);
		return pos_cats;
	}

	/**
	 * Define the negative categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined negative categories
	 */
	protected Map<String, Boolean> defineNegCategories()
	{
		Map<String, Boolean> neg_cats = new HashMap<String, Boolean>();
		neg_cats.put("negative_regulation", false);
		return neg_cats;
	}
	
	/**
	 * Define the neutral categories to be used in the ActivityFlowEdgeOntology
	 * 
	 * @return the defined negative categories
	 */
	protected Map<String, Boolean> defineNeutralCategories()
	{
		Map<String, Boolean> neutral_cats = new HashMap<String, Boolean>();
		neutral_cats.put("regulation", false);
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
			if (afOntology.isDefinedSourceType(conType))
				afCount++;
			if (prOntology.isDefinedSourceType(conType))
				prCount++;
		}
		
		if (afOntology.isDefinedSourceType(refType) && afCount == conCount)
		{
			return afOntology.getDifferentialEdge(refEdge, conEdges, cutoff);
		}
		if (prOntology.isDefinedSourceType(refType) && prCount == conCount)
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
			if (afOntology.isDefinedSourceType(conType))
				afCount++;
			if (prOntology.isDefinedSourceType(conType))
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
	
	@Override
	public boolean isSymmetricalSourceType(String category)
	{
		if (afOntology.isDefinedDiffCategory(category))
			return afOntology.isSymmetricalSourceType(category);
		return prOntology.isSymmetricalSourceType(category);
	}

	@Override
	protected Color getDifferentialEdgeColor(String category)
	{
		if (afOntology.isDefinedDiffCategory(category))
			return afOntology.getDifferentialEdgeColor(category);
		return prOntology.getDifferentialEdgeColor(category);
	}

	@Override
	protected Color getSourceEdgeColor(String edgeType)
	{
		if (afOntology.isDefinedSourceType(edgeType))
			return afOntology.getSourceEdgeColor(edgeType);
		return prOntology.getSourceEdgeColor(edgeType);
	}
	
	@Override
	protected ArrowHead getDifferentialEdgeArrowHead(String category)
	{
		if (afOntology.isDefinedDiffCategory(category))
			return afOntology.getDifferentialEdgeArrowHead(category);
		return prOntology.getDifferentialEdgeArrowHead(category);
	}

	@Override
	protected ArrowHead getSourceEdgeArrowHead(String edgeType)
	{
		if (afOntology.isDefinedSourceType(edgeType))
			return afOntology.getSourceEdgeArrowHead(edgeType);
		return prOntology.getSourceEdgeArrowHead(edgeType);
	}
	
	/**
	 * Provide a default child-parent mapping between source categories.
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
	protected void addDefaultPaintedParents()
	{
		prOntology.addColor("ptm", Color.BLUE);
		prOntology.addArrowHead("ptm", ArrowHead.ARROW);
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

}
